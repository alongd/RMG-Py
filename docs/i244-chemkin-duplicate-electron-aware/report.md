# I-244 — the Chemkin writer marked two reactions of different stoichiometry as duplicates

Worktree `/home/alon/Code/RMG-Py-i244-chemkin-duplicate-electron-aware`, branch
`i244-chemkin-duplicate-electron-aware`, cut from `plasma` at `311818121`.
Database `/home/alon/Code/RMG-database-plasma`, read-only and unmodified.
All runs captured to `logs/` with both streams.

> **Round 2 superseded parts of §§1–8.** A review round found that the refinement below was added
> to one side of a group/pairwise mismatch, and two claims here are false. Read §11 onwards
> alongside them; each superseded statement is flagged in place.

---

## 1. What changed

One function, `mark_duplicate_reaction` in `rmgpy/chemkin.pyx`. The branch that **marks** a pair
now consults the per-side electron placement; the two branches that **un**-mark a pair are
untouched. Nothing else in the tree changed.

> **SUPERSEDED by §11.** Putting the refinement on the marking branch leaves it unreachable for a
> pair that arrives already flagged, which is the state the defect lives in. The repair is
> group-level; `mark_duplicate_reactions` now recomputes every flag from a group key, and
> `save_chemkin` calls it.

```python
electrons1 = get_electron_placement_counts(reaction1)
electrons2 = get_electron_placement_counts(reaction2)
same_dir_duplicate     = same_dir_match     and electrons1 == electrons2
opposite_dir_duplicate = opposite_dir_match and electrons1 == (electrons2[1], electrons2[0])
```

## 2. The rule, and its truth table

> **Two reactions may be marked duplicates only if, in addition to their reactant and product
> lists matching, their per-side electron placement counts match in the same orientation the
> lists matched in** — side for side on a same-direction match, each side against the other
> reaction's opposite side on a reverse match.

That is the same rule `Reaction.is_isomorphic` and `rmgpy.rmg.model.are_identical_species_references`
already apply; this was the one identity comparison on the export path that had not been given it.

| orientation | counts₁ | counts₂ | before | after | reachable as |
|---|---|---|---|---|---|
| same | (0, 0) | (0, 0) | mark | **mark** | all neutral chemistry |
| same | (1, 2) | (1, 2) | mark | **mark** | two copies of a declared ionisation channel |
| same | (0, 1) | (0, 1) | mark | **mark** | two undeclared owners of one channel |
| same | (1, 2) | (0, 1) | mark | **no mark** | declared channel beside an undeclared writing of it |
| same | (1, 0) | (2, 1) | mark | **no mark** | radiative beside three-body recombination |
| opposite | (1, 0) | (0, 1) | mark | **mark** | a genuine forward/reverse pair |
| opposite | (1, 2) | (1, 0) | mark | **no mark** | ionisation beside radiative recombination |

Every row above the reachability column is driven through the real writer and read off the emitted
deck in `test/rmgpy/i244ChemkinDuplicateElectronTest.py::TestPlacementTruthTable` and its siblings.

**It is a strict refinement**, and the shape of the change is what makes it one. The new condition
is a conjunct on the marking branch only, so it can turn a `True` into a `False` and nothing else.
Narrowing `same_dir_match`/`opposite_dir_match` themselves — the obvious-looking alternative — is
*not* a refinement: a pair that arrives already flagged `duplicate=True` and is currently unmarked
by the opposite-direction-irreversible branch would stop matching at all, the branch would never
run, and the pair would keep its `DUPLICATE` lines. That wrong shape is pinned red; see §6.

> **RETRACTED — this paragraph is an arithmetic fact dressed as a safety property.** "Can only turn
> a `True` into a `False`" is true and is not a guarantee of anything: under-marking is the
> direction that makes Chemkin reject a file. See §11.

~~For every reaction outside the plasma families and libraries both counts are `(0, 0)`, so the
marking verdict is unchanged bit for bit.~~

> **RETRACTED — false, measured two ways in §12.** Of 140 families, 17 carry a nonzero `electrons`
> and **11 of them are not plasma families**: six `Cation_*` and five
> `Surface_Proton_Electron_Reduction_*`, all `electrons = -1`, hence counts `(1, 0)`. What is true
> is that all eleven are *one-sided*, so declared and net-derived placement agree and a verdict
> moved for them is moved correctly.

## 3. The reproduction

`docs/i244-chemkin-duplicate-electron-aware/probe_repro.py` — the real loader, a real
`CoreEdgeReactionModel` so `make_new_species` unifies isomorphic species into single objects, then
`render_chemkin_file(..., check_for_duplicates=True)`, which is what `save_chemkin_file` calls.
Unification verified in the log by object identity: both reactions hold the *same* `[Li]` and
`[Lip]` objects.

RMG logs the warning the brief predicted, in the representation that hides the bug:

```
WARNING:root:Marked reaction [Li](2) => [Lip](3) as duplicate of [Li](2) => [Lip](3) for saving to Chemkin file.
```

and the emitted deck (`logs/repro-BEFORE-stdout.log`) contains:

```
REACTIONS    KCAL/MOLE   MOLES

[Li](2)+e-(1)=>[Lip](3)+e-(1)+e-(1)                 1.490e+18 -0.267    162.150
    TDEP/e-(1)/   ! VoronovEIArrhenius exported as ...
DUPLICATE

[Lip](3)+e-(1)=>[Li](2)                             1.734e+14 -0.801    0.010
    TDEP/e-(1)/   ! BadnellRRArrhenius exported as ...

[Li](2)=>[Lip](3)+e-(1)                             1.000000e+10 0.000     0.000
DUPLICATE

END
```

Two reactions stamped `DUPLICATE` with **2 reactants / 3 products** against **1 reactant /
2 products**. Both carry `electrons = +1`; only the per-side pair `(1, 2)` against `(0, 1)`
separates them.

After the fix (`logs/repro-AFTER-stdout.log`): `RESULT: 0 of 3 reactions marked DUPLICATE`, and the
deck carries no `DUPLICATE` line.

Both entries arrive with `duplicate=False` and are overwritten by the writer, which confirms the
brief's third finding — a library cannot defend itself by declaring `duplicate=False`.

## 4. The Cantera leg follows — measured

`rmgpy/yaml_cantera2.py` never calls `mark_duplicate_reactions`; it reads `reaction.duplicate` at
line 681. The probe renders the Cantera YAML from the *same* reaction objects, after the Chemkin
writer has run:

```
==== Cantera YAML reaction entries ====
  equation='[Li](2) + e-(1) => [Lip](3) + e-(1) + e-(1)'   duplicate=False
  equation='[Lip](3) + e-(1) => [Li](2)'                   duplicate=False
  equation='[Li](2) => [Lip](3) + e-(1)'                   duplicate=False
```

Before the fix the first and third carried `duplicate=True`. The Cantera writer needs no change.
Confirmed again end to end: a full `rmg.py` run of the lithium plasma deck
(`docs/i123-integration/input.py`, `generateCanteraYAML2=True`) exits 0 and its
`cantera2/chem_edge_annotated.yaml` contains no `duplicate:` key at all
(`logs/plasma-deck-stdout.log`).

## 5. Non-plasma export is byte-identical

> **SUPERSEDED by §12.** True of the electron refinement alone. It is no longer true of the
> branch: the group recompute deletes exactly one `DUPLICATE` line from the edge deck (11 → 10),
> and that line was a lone `DUPLICATE` that Cantera's parser rejects the base deck for.

`examples/rmg/minimal/input.py` run end to end through `rmg.py`, once before the change and once
after, into separate directories:

```
$ diff -r <before>/chemkin <after>/chemkin
$ echo $?
0
```

Zero lines of difference across every deck in the directory — core and edge, every iteration
snapshot. The edge deck carries **11 `DUPLICATE` lines before and 11 after**, so real non-plasma
duplicate marking is demonstrably still happening, not merely absent from both sides.
(`logs/nonplasma-chemkin-diff.txt`, empty.)

The Cantera directory differs only in the generation timestamp
(`logs/nonplasma-cantera-diff.txt`):

```
4c4
< date: Mon, 21 Sep 2026 10:03:50 +0300
---
> date: Mon, 21 Sep 2026 10:05:46 +0300
```

## 6. Tests, and every red state

`test/rmgpy/i244ChemkinDuplicateElectronTest.py`, 13 cases. Every assertion is on the text of the
emitted deck, parsed out of the REACTIONS block, except one that is explicitly labelled — see §7.

Four separate inversions, each built and run, together drive **all 13** red:

| inversion | logged | tests it reds |
|---|---|---|
| fix removed (`same_dir_duplicate = same_dir_match`) | `logs/RED1-fix-removed-*.log` | 4 — the two placement-mismatch cases, the three-body separation, both `False` truth-table rows |
| **wrong shape** — narrow the match flags themselves | `logs/newtests-RED2-wrongshape-*.log` | 2 — both un-marking branches |
| over-refined (`same_dir_duplicate = False`) | `logs/RED3-overrefined-*.log` | 6 — every "still a duplicate" case, incl. all three `True` truth-table rows |
| orientation guard removed from `expand_electrons` | `logs/RED4-orientation-guard-removed-*.log` | 1 — the export-boundary pin |

The second row is the one worth reading: those two tests exist *because* the plausible fix is wrong,
and they catch it.

Green after restoring, with the `.so` rebuilt: **13 passed** (`logs/newtests-green2-stdout.log`).

## 7. Suites — my own measurements

Interpreter pinned to `/home/alon/anaconda3/envs/rmg_env/bin/pytest`.

| suite | result |
|---|---|
| `test/rmgpy/chemkinTest.py` + `reactionTest.py` + `yaml_cantera1Test.py` + `yaml_cantera2Test.py` + the new file | **218 passed, 12 skipped** |
| `test/rmgpy/rmg/` | **149 passed, 2 skipped, 1 failed + 1 error** (see below) |
| all of `test/rmgpy/*.py` and `test/rmgpy/rmg/`, `-m "not functional and not database"` | **1106 passed, 4 skipped, 102 deselected, 1 failed + 1 error** |

The single failure is `mainTest.py::TestProfiling::test_make_profile_graph`, and it is
environmental and pre-existing:

```
subprocess.CalledProcessError: Command '['ps2pdf', '.../RMG.profile.dot.ps2', '.../RMG.profile.dot.pdf']'
    returned non-zero exit status 1
Error: /undefinedfilename in (.../RMG.profile.dot.ps2)
```

`which dot` finds nothing — graphviz is not installed in this environment — so the `.ps2` is never
produced and `ps2pdf` is handed a missing file. `ps2pdf` itself works (`ps2pdf t.ps2 t.pdf` exits 0
on a trivial file). The test calls `make_profile_graph` on a committed `.profile` and touches no
reaction model, no writer, and no Chemkin code.

## 8. Contradictions with the brief — the measurement won three times

**(a) The suggested reproduction does not reproduce.** The brief says to build the reproduction
from `PlasmaRadiativeRecombination` `(1, 0)` and `PlasmaElectronImpactIonization` `(1, 2)`. Loaded
through the real path with unification, those two shipped entries are `[Lip] => [Li]` and
`[Li] => [Lip]`, and **both declare `reversible = False`**. That is an opposite-direction,
both-irreversible pair, which `mark_duplicate_reaction` already declines to mark — for a reason
that has nothing to do with electrons. Measured: `RESULT: 0 of 2 reactions marked DUPLICATE`.

The defect needs a **same-direction** pair. Structurally it has to: if two reactions have the same
heavy species on each side and both balance in the `E` pseudo-element, their *net* electron counts
are necessarily equal, so a same-direction false duplicate can only arise where a **two-sided
declaration** makes the incident counts differ at equal net. On this branch that is
`PlasmaElectronImpactIonization`'s `(1, 2)` beside the net-derived `(0, 1)` — which is exactly the
shape I-148 recorded from a seed round trip ("the ionisation channel's `(1, 2)` collapsed to
`(0, 1)`, the restarted core carried the channel twice") and exactly what any user library that
writes the channel without declaring its incident order produces. That is the reproduction I built.

**(b) `rmgpy/rmg/yaml_cantera2.py` does not exist.** The Cantera writer is at
`rmgpy/yaml_cantera2.py`. The line number (681, `if reaction.duplicate:`) is right.

**(c) The three findings.** Findings 2 and 3 confirmed as stated. Finding 1 — "the net electron
count is `-1` for *both* reactions" — is right in substance and wrong in the specific number for
the case reachable here: in my reproduction the net is `+1` for both. The load-bearing part of the
finding, that the net scalar is *equal* for the two and therefore cannot separate them, holds
exactly.

**(d) A fourth finding, not in the brief.** The two halves of `rmgpy/electron_balance.py` disagree
about reversed orientation. `get_electron_placement_counts` reads a declaration in either
orientation (returning `(2, 1)` for a reaction stored backwards against a `(1, 2)` declaration);
`expand_electrons` accepts only the forward one and raises `MechanismWriterError` otherwise. The
asymmetry is deliberate and documented on the export side, but it has a consequence worth naming:
a three-body recombination **cannot** reach a deck by borrowing the ionisation family's declaration
backwards. Pinned in `test_the_export_boundary_refuses_a_reversed_declared_reaction`.

## 9. The argon three-body branches — is this fix sufficient?

The brief asks explicitly. The answer is **yes at the identity level, conditional on one thing**,
and it will need re-verification when those branches land.

- The pair the brief reproduced — `[Arp] + e- => [Ar]` placement `(1, 0)` against
  `[Arp] + 2e- => [Ar] + e-` placement `(2, 1)` — is a same-direction match with differing
  placements, so the refined rule declines to mark it. Measured on the lithium analogue in
  `test_the_duplicate_check_separates_the_pair`.
- **The condition:** the three-body library must carry its own **forward** `(2, 1)` entry in
  `FAMILY_ELECTRON_PLACEMENT`. If it does not, its reactions fall back on the net-derived rule,
  `electrons = -1` becomes `(1, 0)`, and they are indistinguishable from radiative recombination —
  still marked duplicate, and, more seriously, still *exported* second order with the incident
  electron lost. That limit is pinned rather than left to be discovered, in
  `test_an_undeclared_three_body_owner_is_still_folded_onto_the_radiative_channel`.
- Per §8(d), the three-body owner also cannot borrow `(1, 2)` reversed: `expand_electrons` refuses
  it at the export boundary. So the declaration has to be its own.

**Recommendation:** when the two unmerged branches land, re-run this ticket's tests against them
and check that the three-body library's owner is declared `(2, 1)` forward. I did not open, modify,
or read into those branches.

## 10. Housekeeping

- No rate, library entry, training entry, family, or database file was modified. The
  RMG-database checkout is untouched.
- `rmgrc` was created from `rmgrc.template` (git-ignored, as designed).
- `git status` is clean apart from the intended change and the new files.
- Nothing pushed, merged, or opened as a PR.

---

# Round 2 — `duplicate` is a group flag, and the refinement sat on one side of that

## 11. Defect 1 — reproduced both ways before anything was changed

Two engines, both already built, neither one rebuilt or reverted to get the other:
the branch head `93e070b0f` in this worktree, and its base `311818121` at
`/home/alon/Code/RMG-Py-plasma`. `docs/i244-chemkin-duplicate-electron-aware/defect1_repro.py`
takes the worktree as `argv[1]` and `chdir`s into it **before** importing `rmgpy`, because the
shared editable install resolves `import rmgpy` from the current working directory. Each run
prints the loaded `.so` path and a value-level fingerprint (the branch's
`mark_duplicate_reaction.__doc__` contains "Electrons are participants"; the base's does not), so
the two arms cannot be confused. Logs: `logs/D1-repro-BASE-*`, `logs/D1-repro-BRANCH-*`.

### 11a. The pre-marked bypass

The refinement lives in the `else` of `if reaction1.duplicate and reaction2.duplicate`, so a pair
that arrives already flagged never reaches it. Driving `mark_duplicate_reaction` with the exact
pair this ticket exists to separate — declared `(1, 2)` against net-derived `(0, 1)`, both
`electrons = +1`:

| arriving flags | counts | base `311818121` | branch `93e070b0f` |
|---|---|---|---|
| `(True , True )` | `(1,2)` vs `(0,1)` | `(True, True)` | **`(True, True)`** |
| `(True , False)` | `(1,2)` vs `(0,1)` | `(True, True)` | `(True, False)` |
| `(False, True )` | `(1,2)` vs `(0,1)` | `(True, True)` | `(False, True)` |
| `(False, False)` | `(1,2)` vs `(0,1)` | `(True, True)` | `(False, False)` |

The first row is the defect: the same two reactions, the same placements, shipped as duplicates of
each other. Taken all the way to a deck, the base emits both stamped `DUPLICATE` for **every**
arriving combination, and so does the branch for the pre-marked one:

```
BASE  arrive (True, True)   -> Li+e-=>Lip+e-+e- DUPLICATE=True, Li=>Lip+e- DUPLICATE=True
BASE  arrive (False, False) -> Li+e-=>Lip+e-+e- DUPLICATE=True, Li=>Lip+e- DUPLICATE=True
```

### 11b. The cross-group clearing

Four legitimately pre-marked reactions in two groups — two `(1, 2)` not pressure dependent, two
`(0, 1)` pressure dependent. A cross-group comparison reaches the mixed-pressure-dependence
un-marking branch, which clears **both** flags, one of which belonged to the group that was never
in question:

```
arriving   A1=True A2=True B1=True B2=True     placements A1,A2=(1,2)  B1,B2=(0,1)
resulting  A1=False A2=False B1=True B2=True
deck       Li+e-=>Lip+e-+e-  DUPLICATE=False   <-- written twice, neither marked
           Li+e-=>Lip+e-+e-  DUPLICATE=False       CHEMKIN REJECTS THIS DECK
           Li(+M)=>Lip+e-(+M) DUPLICATE=True
           Li(+M)=>Lip+e-(+M) DUPLICATE=True
```

**Correction to the brief.** This half is *identical on the base engine* — same four flags, same
deck. It is a pre-existing defect that the electron refinement neither introduced nor worsened;
the brief presents it as the converse of the refinement, which it is not. My own first trace
predicted the two arms would differ, and the measurement refuted that too: the marking branch
already tests pressure dependence, so the cross-comparison that I expected to re-mark `A1` never
fires.

## 12. The ending I chose, and what it costs

> **Group-level recomputation.** `mark_duplicate_reactions` no longer sweeps pairs. It partitions
> the list by a new `chemkin_duplicate_group_key` — class, collider identity, the ordered
> participant identities of each side **paired with that side's free-electron count**, pressure
> dependence, and reversibility (direction-specific when irreversible, the unordered pair of sides
> when reversible) — and sets `duplicate = len(group) > 1` for every reaction. `save_chemkin` calls
> it once over exactly the list it is about to write.

Why not the other two endings the brief offers. **Refining the un-marking branches** cannot work:
it makes those branches stop firing for a placement-mismatched pre-marked pair, which leaves the
pair's flags exactly as they arrived — `(True, True)` — and two lone `DUPLICATE` lines in the deck.
It moves the failure, it does not remove it. **Refusing** on a pre-marked pair whose placements
differ would reject §11b's four-reaction model, which is *legal input*: two groups whose members
legitimately have different placements. The rule would fire on correct mechanisms.

What group recomputation costs, in order of how much it should worry a reviewer:

1. **It can clear a flag, and that reaches non-plasma decks.** A reaction that arrives flagged and
   has no group mate in the list being written now loses its flag. Nothing before ever did that.
   Measured on `examples/rmg/minimal`, base against fixed, same input and same database: the
   chemkin diff is **one deleted `DUPLICATE` line** per affected edge deck and nothing else — no
   reordering, no rate change, no lost pair. The line removed was
   `C4H7(72)<=>C4H7(44)`, which appears **once** in a 398-reaction deck.

   Cantera's own Chemkin parser settles what that line was worth. On the base deck:

   ```
   InputFileError thrown by Kinetics::checkDuplicates:
   No duplicate found for declared duplicate reaction number 179 (C4H7(72) <=> C4H7(44))
   ```

   On the fixed deck, `Validating mechanism... PASSED`. So RMG has been emitting an **invalid edge
   deck on plain neutral chemistry**, and this repair is what makes it valid. That is a wider
   finding than the ticket and it is the strongest evidence the clearing is right rather than
   risky. Logs: `logs/ck2yaml-BASE-stdout.log`, `logs/ck2yaml-FIXED-stdout.log`,
   `logs/nonplasma-chemkin-diff-groupfix.txt`.

2. **It changes `save_chemkin`,** which is the one scope widening I took. Every `render_*` call
   inside `save_chemkin` passes `check_for_duplicates=False` with the comment *"We should already
   have marked everything as duplicates by now"* — so without a call there the repair would never
   reach a real RMG run's deck, only direct `save_chemkin_file` callers. One line, and the comment
   now states what is true instead of what was hoped.

3. **The pairwise `mark_duplicate_reaction` is left as it is**, deliberately, and is now documented
   as an incremental hint rather than an authority. `rmgpy/rmg/model.py` calls it while enlarging
   the model, where `reaction_list` is only what has been checked so far; clearing a flag on that
   evidence would be wrong, and would feed `check_for_existing_reaction`. Its marking condition now
   asks `chemkin_duplicate_group_key` the same question the recompute asks, so the two cannot drift
   about what "the same entry" means, but it can still under-mark and its docstring says so
   plainly.

4. **Cost in time is negative.** One hashing pass, linear, against the previous quadratic pairwise
   sweep. The brief anticipated a cost here; there is not one.

The docstring's old safety argument — "a conjunct on the marking branch can only turn a `True` into
a `False`, so the change can only mark fewer pairs" — is rewritten. It was arithmetic, not a safety
property: under-marking is the direction Chemkin rejects. What the change now claims is only that
the pairwise function never introduces a mark the group recompute would not also make, and that
this is bookkeeping rather than correctness.

## 13. Defect 2 — the universal negative, measured two ways

Both methods run by me against `/home/alon/Code/RMG-database-plasma`, and they agree.
Logs: `logs/D2-census-static-*`, `logs/D2-census-dynamic-*`.

| | families | nonzero `electrons` | plasma-named | **not** plasma-named |
|---|---|---|---|---|
| static grep over `*/groups.py` | 140 | 17 | 6 | **11** |
| dynamic census over loaded families | 140 | 17 | 6 | **11** |

The eleven: `Cation_Addition_MultipleBond`, `Cation_Addition_MultipleBond_Disprop`,
`Cation_Li_Abstraction`, `Cation_NO_Ring_Opening`, `Cation_NO_Substitution`,
`Cation_R_Recombination`, `Surface_Proton_Electron_Reduction_Alpha`,
`…_Alpha_vdW`, `…_Beta`, `…_Beta_Dissociation`, `…_Beta_vdW` — every one `electrons = -1`,
counts `(1, 0)`.

So the claim "for every reaction outside the plasma families and libraries the counts are `(0, 0)`
on both sides" is **false**. The corrected statement, which is narrower and is what actually bounds
the change: 123 of 140 families key `(0, 0)` and cannot move; the other 17 place **one-sidedly**
except one, so declared and net-derived placement agree for them and a verdict moved is moved
correctly; and only a **two-sided** declaration can produce an answer the net count could not.
There are exactly two in the database — the family `Plasma_Electron_Impact_Ionization` and the
library `PlasmaElectronImpactIonization`, both `(1, 2)`.

The false gloss had spread to three more places, all corrected in this commit rather than only at
the location the round named: `rmgpy/chemkin.pyx` (`chemkin_duplicate_group_key`),
`rmgpy/rmg/model.py::are_identical_species_references`, and
`test/rmgpy/i134DuplicateElectronsTest.py::TestTheUndeclaredOwnerContract`. The commit body of
`5db81eef0` carries it too and cannot be corrected without rewriting that commit; it is retracted
here and in this report instead.

## 14. Defect 3 — the claim is refuted; the repoint was still necessary

The brief says the test at `i244ChemkinDuplicateElectronTest.py:438` "uses mismatched reactions
that have no valid duplicate mate, so the branch it names cannot fire and the test passes against
the broken behaviour."

**Measured, the branch fires.** Running that test with live log capture against the unmodified
branch build prints the branch's own warning, which no other code emits
(`logs/D3-line438-warnings-stdout.log`):

```
WARNING root: Marked reaction Li(2) => Lip(3) as not duplicate because of mixed pressure
              dependence for saving to Chemkin file.
```

And building the pre-fix code with both un-marking branches disabled makes the original,
deck-level form of *both* tests fail (`logs/D3-oldform-prefixcode-branch-disabled.log`):

```
OLD-FORM opposite_direction_irreversible_branch_still_unmarks : FAIL
OLD-FORM mixed_pressure_dependence_branch_still_unmarks       : FAIL
```

So against the code they were written for, they could fail and they were not protecting the broken
behaviour.

**The real finding is the converse, and it arrives with this commit.** Under the group recompute,
the same original deck-level assertions pass with the branches deleted
(`logs/D3-oldform-newcode-branch-disabled.log`):

```
OLD-FORM opposite_direction_irreversible_branch_still_unmarks : PASS
OLD-FORM mixed_pressure_dependence_branch_still_unmarks       : PASS
```

— because the recompute produces the same deck either way. Both tests would have become
ninth-instance checks that cannot fail, created by this very repair. Both are therefore repointed
to drive `mark_duplicate_reaction` directly and assert the branch-specific warning text, and both
are shown red against a build with those branches disabled
(`logs/RED-unmarking-removed-stdout.log`, 2 failed / 17 passed).

## 15. Which tests are load-bearing

Removing the electron term from the group key (`(0, 0)` forced) and rebuilding
(`logs/RED-refinement-removed-stdout.log`) turns **7 of 19** red:

1. `TestElectronPlacementSeparatesDuplicates::test_two_owners_of_one_ionisation_channel_are_not_duplicates`
2. `TestThreeBodyRecombinationShape::test_the_duplicate_check_separates_the_pair`
3. `TestPlacementTruthTable::test_same_direction_rows[(1, 2) vs (0, 1)]`
4. `TestPlacementTruthTable::test_same_direction_rows[(0, 1) vs (1, 2)]`
5. `TestDuplicateIsAGroupPredicate::test_a_pre_marked_placement_mismatched_pair_is_cleared`
6. `TestDuplicateIsAGroupPredicate::test_the_verdict_does_not_depend_on_the_flags_it_arrives_with`
7. `TestDuplicateIsAGroupPredicate::test_the_production_save_path_recomputes_the_flags`

Items 1–4 are exactly the four the brief predicted among the original thirteen — that figure is
confirmed. The other nine originals are controls or pin pre-existing behaviour, and three of my six
new tests are load-bearing for the group repair rather than for the electron term
(`…does_not_clear_another_groups_flags`, `…lone_flag_that_names_no_mate_is_cleared`,
`…mixed_reversibility_pair_is_cleared`); they go red against the pre-fix build, which is where
their red state was demonstrated.

Every new test sets `duplicate = True` on at least one input before calling. All six were shown red
against the pre-fix build (`logs/D1-newtests-RED-stdout.log`, 6 failed / 13 passed) and green after
(`logs/D1-newtests-GREEN-stdout.log`, 19 passed).

## 16. Suites

| suite | result |
|---|---|
| `test/rmgpy/i244ChemkinDuplicateElectronTest.py` | 19 passed |
| `chemkinTest` + `plasmaExportTest` + `i134` + `i126` + `i135` + `reactionTest` + `i209` | 457 passed, **5 failed**, 2 skipped, 1 xfailed |
| `test/rmgpy` (unit, excluding `test/rmgpy/rmg`) | 2963 passed, 32 skipped, 121 deselected |
| `test/rmgpy/rmg` | 150 passed, 2 skipped |

The 5 failures are all in `i134DuplicateElectronsTest.py` and are **pre-existing**: run against the
base engine at `/home/alon/Code/RMG-Py-plasma` they fail identically, 5 failed / 178 passed
(`logs/REG-i134-BASE-stdout.log`). They are database-content failures —
`PlasmaRadiativeRecombination has 2 entries, expected 1`, the second being the argon entry a later
branch added — and have nothing to do with duplicate marking.

## 17. Housekeeping, and one false positive recorded rather than worked around

- `git diff --check` over the already-committed logs went from **916** flagged lines to **0**. The
  916 were trailing whitespace inside committed evidence logs. Rewriting a log to satisfy a
  whitespace check would falsify the evidence, so a new `.gitattributes` marks `docs/**/logs/**`
  `-whitespace` instead.
- **Two** flagged lines survive, both in logs added by this commit, and both are false positives in
  the check itself: pytest's `======= 5 failed, 457 passed … ========` summary banner is read as a
  *leftover conflict marker*. Suppressing that needs `-diff`, which would make every evidence log
  show as "binary" in review. Recorded here as a finding rather than worked around.
- A previously-committed log, `logs/FINAL-verifier-stdout.log`, was overwritten by a round-2 run
  before I noticed that `tee` truncates. It has been restored from `93e070b0f`; this round's run is
  `logs/FINAL-verifier-round2-stdout.log`.
- Two stray files, `chem-gas.yaml` and `chem_annotated.yaml`, appeared at the repo root during the
  test runs (a suite writes them to the cwd) and were deleted. They were created by this session,
  not pre-existing.
- `docs/pm-verify/` and every other untracked path found at the start is preserved.
- No rate, library entry, training entry, family, or database file touched. Nothing pushed,
  merged, or opened as a PR.
