# I-244 round 75 — three HIGH, one MEDIUM, one LOW

Measured on `4be053318`, branch `i244-chemkin-duplicate-electron-aware`, nothing pushed.

Every claim below ends at `cantera.Solution`. Deck text and `ck2yaml.convert_mech` are both
incapable of seeing this defect class: `convert_mech` writes YAML without running
`Kinetics::checkDuplicates`, so the measured signature is **`convert_mech` ACCEPTED;
`Solution` REJECTED**, and three of the five findings produce exactly that.

---

## 1. The enumeration, put first because it is what the round asked for

The interface the defect travels through is *"anything that emits a Chemkin `DUPLICATE`
line or a Cantera `duplicate: true` key"*. Rounds 73 and 74 repaired against a list of that
interface that was incomplete, and verified against a subset of the incomplete list. Here is
the whole census, built by grepping every `.duplicate` access and every registered writer
rather than by recalling which ones mattered.

### Registered writers (`rmgpy/rmg/main.py`)

| # | Writer | Registered | Emits duplicates? | State after round 75 |
|---|---|---|---|---|
| 1 | `ChemkinWriter` | `main.py:904` | yes | recomputes per list (round 73) |
| 2 | `RMSWriter` | `main.py:911` | **no** — RMS has no duplicate concept | N/A, confirmed by grep |
| 3 | `CanteraWriter1` | `main.py:913` | yes | **repaired here** — was reading the flag |
| 4 | `CanteraWriter2` | `main.py:919` | yes | recomputes per list (round 74) |
| 5 | `OutputHTMLWriter` | `main.py:925` | no mechanism artifact | N/A |
| 6–8 | `QMDatabaseWriter`, `ExecutionStatsWriter`, `SimulationProfileWriter` | — | no | N/A |

### Every site that decides a duplicate answer

| # | Site | Reached from | Round-75 state |
|---|---|---|---|
| 1 | `chemkin.save_chemkin_file` | ChemkinWriter, `scripts/` | recomputes |
| 2 | `chemkin.save_chemkin_surface_file` | ChemkinWriter | recomputes |
| 3 | `chemkin.mark_duplicate_reactions` | `model.py:2182` | recomputes, mutates deliberately |
| 4 | `chemkin.write_kinetics_entry` (MULTI leaves) | every Chemkin path | **HIGH 3 — repaired** |
| 5 | `chemkin.mark_duplicate_reaction` | `model.py:833`, `:1962` | **HIGH 2 — repaired** |
| 6 | `yaml_cantera2._collect_reaction_entries` | CanteraWriter2 | recomputes |
| 7 | `yaml_cantera2.reaction_to_dict_list` (MULTI leaves) | CanteraWriter2 | **HIGH 3 — repaired** |
| 8 | `yaml_cantera2.get_reaction_equation` (collider) | CanteraWriter2 | **MEDIUM — repaired** |
| 9 | `yaml_cantera1._collect_reactions` | CanteraWriter1 | **HIGH 1 — repaired** |
| 10 | `yaml_cantera1._build_equation_string` (collider) | CanteraWriter1 | **MEDIUM — repaired (not in the round's list)** |
| 11 | `reaction.Reaction.to_cantera` (MULTI leaves) | Writer1, `tools/canteramodel.py` | **HIGH 3 — repaired (not in the round's list)** |
| 12 | `reaction.Reaction.to_cantera` (single) | Writer1, `tools/canteramodel.py` | reads `self.duplicate`; closed for Writer1 by #9, open for `canteramodel` — DEFERRED |
| 13 | `arkane/pdep.py:506` | Arkane | one reaction, no list to key over — DEFERRED, own ticket |
| 14 | `rmgpy/rmg/output.py:647` | `diffModels.py` HTML | same false premise as HIGH 3, but display only — enumerated, not repaired |
| 15 | `data/kinetics/library.py`, `common.py` | library load | input side, not a writer |

**Two of those rows were not in the round's findings and are load-bearing.** Row 11 is the
site `CanteraWriter1` and `rmgpy/tools/canteramodel.py` actually go through, so repairing
only the two named HIGH-3 sites would have left Writer1 emitting the same invalid mechanism.
Row 10 is the same MEDIUM in the other writer.

---

## 2. HIGH 1 — the leak is relocated, not closed. Confirmed and closed.

**Reproduced** (`logs/R75-RED-findings-stdout.log`):

```
production flags over core+edge: [True, True]
group answer for the CORE-ONLY list: [False]
FAIL  Writer1, core-only export after production marking
      InputFileError thrown by Kinetics::checkDuplicates:
      No duplicate found for declared duplicate reaction number 0 (H(1) + OH(2) => H2(3) + O(4))
PASS  Writer2, core-only export after production marking
```

`yaml_cantera1._collect_reactions` now keys its own list with `chemkin_duplicate_flags`,
exactly as `yaml_cantera2._collect_reaction_entries` does, and `reaction_to_dicts` writes
that answer into the entry rather than letting `Reaction.to_cantera`'s value stand.

The round asked for the alternative — evidence that Writer1 cannot receive a
production-marked reaction — and it does not exist: `CanteraWriter1.update` is handed
`rmg.reaction_model.core.reactions`, the same shared objects `model.py:833` marks over
core+edge.

## 3. HIGH 2 — cross-class through the pairwise path. Confirmed and closed, with the decision stated.

**Reproduced:**

```
production (pairwise) flags: [False, False]
group answer over the same list: [True, True]
FAIL  Writer1, cross-class pair
      Undeclared duplicate reactions detected
```

**The decision the round asked me to make explicitly: the class comparison goes, AND the
writers must not trust pairwise output. Both, because neither alone is sufficient.**

- Removing the term alone still leaves a writer serialising a hint computed over some other
  list — that is HIGH 1, and it is a different defect.
- Making every writer recompute alone still leaves `Reaction.duplicate` wrong for everything
  that reads it without recomputing: `arkane/pdep.py`, `tools/canteramodel.py`, and the HTML
  report.

The term's provenance argues for keeping it, and is why the round asked. It is stock RMG's,
and it carries its author's own unanswered doubt inline: *"RHW question: why can't
TemplateReaction be duplicate of LibraryReaction, in Chemkin terms? I guess it shouldn't
happen in RMG."* The answer is that in Chemkin terms they certainly can — a class is not
serialized, so two reactions of different classes over the same participants write one and
the same equation. Round 74 removed the term from the group key for that reason; leaving it
in the pairwise marker made the two functions disagree about what "the same entry" means,
and model growth runs through the pairwise one. The comment is preserved verbatim in the
docstring rather than deleted with the code.

Note the asymmetry that made this safe to remove: in the pairwise predicate a class mismatch
means `continue` — never marked — while in a group key it means *cleared*. Removing it can
only merge groups, never split them, and under-marking is the direction Chemkin and Cantera
reject.

## 4. HIGH 3 — single-entry `Multi` wrappers. Confirmed at THREE sites, not two.

**Reproduced for both constructors, through all three writers:**

```
-- MultiArrhenius([one]) -- leaves: 1
FAIL  Chemkin -> ck2yaml -> Solution   No duplicate found for declared duplicate reaction number 0
FAIL  Writer2                          (same)
FAIL  Writer1                          (same)
-- MultiPDepArrhenius([one]) -- same three
-- control: MultiArrhenius([two]) -- all three PASS
```

The repair is **not** "a one-leaf wrapper is never a duplicate" — that would swing the
defect to the other side. It is "a one-leaf wrapper falls back to the group answer", so a
one-leaf wrapper that genuinely shares its equation with another reaction is still marked.
Both controls are pinned as tests.

The third site, `rmgpy/reaction.py:565`, is the one the round's enumeration missed, and it
is the one Writer1 and `rmgpy/tools/canteramodel.py` reach.

## 5. MEDIUM — unsupported specific colliders. Confirmed, and the fix is narrower than the round proposed.

**Reproduced:** two PLOG reactions differing only by `specific_collider` key apart
(`[False, False]`), serialize to the identical equation, and Cantera rejects with
*"Undeclared duplicate reactions detected"*.

The round's proposed fix — *"Writer2 should refuse, as Chemkin correctly refuses first"* —
would have removed a working capability, and measuring Chemkin rather than assuming it is
what caught that:

```
Arrhenius            Chemkin REFUSED ChemkinError
MultiArrhenius(2)    Chemkin REFUSED ChemkinError
PDepArrhenius        Chemkin REFUSED ChemkinError
ThirdBody            Chemkin REFUSED ChemkinError      <-- but Writer2 renders this correctly
Lindemann            Chemkin WROTE   H(1)+OH(2)(+Ar(5))=>H2(3)+O(4)(+Ar(5))
```

Chemkin refuses a specific collider on **everything except Troe/Lindemann, `ThirdBody`
included**. Writer2 renders `ThirdBody`+collider as `H + H + Ar => H2 + Ar` — the collider
explicit on both sides — which Cantera loads, and which Chemkin's `(+species)` syntax cannot
express at all. Matching Chemkin exactly would have deleted that.

So the refusal is exactly the set the equation builder silently drops: a collider on
kinetics that is not `ThirdBody`, `Lindemann` or `Troe`. Applied to **both** Cantera
writers. A control test pins that two `ThirdBody` reactions with different colliders still
render distinctly and load.

### What this narrows, stated rather than hidden

A `specific_collider` can no longer ride a `Multi*` wrapper through Writer2, because the
wrapper expands into `Arrhenius`/`PDepArrhenius` leaves. That made
`plasmaExportTest.py::test_every_component_keeps_index_collider_pairs_and_electrons` fail,
and the fixture turned out to be **pinning the defect as expected behaviour**: it asserted
the collider came back in the entry's *note*, i.e. that the equation had dropped it. The
same shape is refused outright by the Chemkin writer, so no RMG run could ever have carried
such a reaction to a deck. Same class as round 74's invalid golden — a reference captured
from the implementation inherits the implementation's bug and then outranks the correction.
The fixture's live coverage (index, flux pairs, electrons surviving reconstruction) is kept;
the collider assertion is replaced by a test that the refusal fires and names
`MultiArrhenius` — the type the caller wrote — rather than the `Arrhenius` of a leaf nobody
constructed.

## 6. LOW — the tie-back was source-text theatre. Replaced.

`test_the_production_marking_this_rests_on_still_exists` asserted
`source.count("mark_duplicate_reaction(rxn, checked_reactions)") == 3`. The round is right
that dead code and a comment both satisfy it; it is also red on a pure rename that changes
nothing. Replaced by two tests that run production:

- `test_the_production_marking_this_rests_on_really_marks` calls
  `CoreEdgeReactionModel.mark_chemkin_duplicates()` — real production code, reachable
  without a database — and asserts the flags it leaves on the shared objects. `_production_mark`
  now calls that too, instead of a hand-written copy of the loop in `enlarge()`.
- `test_the_pairwise_marking_production_uses_agrees_with_the_group_key` pins the function
  the two remaining pairwise call sites call. Both sites sit inside database-dependent
  methods (`enlarge`, `add_reaction_library_to_edge`), so a unit test cannot reach the call;
  it pins the callee's behaviour instead of the caller's text.

The precondition assertion inside the scope test is kept, as asked.

---

## 7. RED arms — every repair broken one at a time

`docs/i244-chemkin-duplicate-electron-aware/round75_red_arms.py`, log
`logs/R75-red-arms-stdout.log`. Each arm mutates one file, rebuilds where the file is
compiled, runs the guarding suite, and asserts a **named** test goes red — `must_fail`, so
an inert mutation cannot be scored as caught.

| Arm | Breaks | Guard that went red | Discriminator flipped |
|---|---|---|---|
| D5 | Writer1 reads the flag again | `test_writer1_core_only_export_is_loadable_after_production_marking` (+1 more) | `writer1_keys_its_own_list` |
| D6 | class comparison back in the pairwise marker | `test_the_pairwise_marking_production_uses_agrees_with_the_group_key` | `pairwise_marker_ignores_class` |
| D7 | Chemkin marks every leaf | `test_the_chemkin_deck_cantera_converts_also_loads` ×2 | `chemkin_one_leaf_writes_no_duplicate` |
| D8 | Writer2 marks every leaf | `test_both_cantera_writers_load` ×2 | `writer2_one_leaf_not_marked` |
| D9 | `to_cantera` marks every leaf | `test_to_cantera_alone_does_not_mark_a_one_leaf_wrapper` | `to_cantera_one_leaf_not_marked` |
| D10 | Writer2 drops the collider | `test_a_grouped_reaction_with_an_unrenderable_collider_is_refused` | `writer2_refusal_names_the_wrapper` |
| D11 | Writer1 drops the collider | `test_writer1_refuses_a_collider_it_cannot_put_in_the_equation` | `writer1_refuses_unrenderable_collider` |

Restored: nine discriminators all true, 110 passed.

### Two traps this harness caught, both in itself

1. **D9 is inert through every writer.** `yaml_cantera1.reaction_to_dicts` now overwrites
   whatever `to_cantera` set, so breaking `reaction.py` changes no writer's output. Without
   `must_fail` naming a guard that reaches `to_cantera` directly, that repair would have been
   scored as tested when nothing tested it. `test_to_cantera_alone_does_not_mark_a_one_leaf_wrapper`
   exists for that reason and asserts on `to_cantera`, not on a writer downstream of it.
2. **The first run of these arms printed `<probe failed>` for every discriminator.** The
   probe's fixture was a unimolecular reaction carrying a bimolecular `A`, which trips a
   units assertion inside `write_kinetics_entry` before anything is rendered. All seven arms
   had already "passed" with no build check at all — the safeguard against a build that did
   not take was itself dead, and silently. Fixed and re-run; the `<probe failed>` run is not
   the evidence.

---

## 8. Verification

### End to end, `examples/rmg/minimal`, with BOTH Cantera writers enabled

Round 74's end-to-end run exercised one writer. This one sets
`generateCanteraYAML1=True` and `generateCanteraYAML2=True`, and then **loads every
mechanism the run produced** — not the final one, every iteration's, core and edge
(`logs/R75-minimal-load-all-stdout.log`):

```
cantera1: 48 files, 76 'duplicate: true' lines total
cantera2: 48 files, 76 'duplicate: true' lines total
chem_annotated.inp:      DUPLICATE=0  -> LOADED (26 species, 66 reactions)
chem_edge_annotated.inp: DUPLICATE=10 -> LOADED (172 species, 398 reactions)

LOADED 98, REJECTED 0
```

The two writers independently produce **the same** 76 marks over the same 48 mechanisms.
That equality is the cross-check that matters: before this round Writer1's answer came from
a different source entirely (the object's flag), and there was nothing anywhere comparing
them. The Chemkin figures are unchanged from round 73 — edge `DUPLICATE=10` and loading,
against base `DUPLICATE=11` and rejected.

### Suites

| Suite | Result |
|---|---|
| `test/rmgpy/i244ChemkinDuplicateElectronTest.py` | 40 passed |
| `test/rmgpy/*.py` (top level) | **5 failed, 1082 passed** — all five in `i134DuplicateElectronsTest.py` |
| `test/rmgpy/rmg/` + `test/arkane/` | 381 passed, 9 skipped |
| `test/rmgpy/data/` | 389 passed, 8 skipped |
| RED arms | 7 arms, every one takes its named guard red; restored 110 passed, 9/9 discriminators true |

The five failures are the known I-255 five, identical on the base commit, owned by another
worktree, and deliberately not touched from here.

**One harness trap worth recording**, because it produced 25 phantom failures: running
`test/arkane/` without `/home/alon/anaconda3/envs/rmg_env/bin` on `PATH` fails 25 tests with
`[Errno 2] No such file or directory: 'symmetry'`. That is a missing external binary, not a
regression — the same shape as round 74's transport-data trap, a rejection whose stated cause
is not its real one. With the env's `bin` on `PATH`: 381 passed, zero failed.

## 9. Still open, deliberately

- **`arkane/pdep.py:506`** — writes one network reaction at a time with no list to key over.
  Own ticket. `arkane/kinetics.py` is NOT a `DUPLICATE` writer; that half of an earlier
  relayed claim stays refuted.
- **`rmgpy/reaction.py:570`** — `to_cantera` still reads `self.duplicate` for a single
  reaction. Closed for `CanteraWriter1` by the Writer1 repair, still open for
  `rmgpy/tools/canteramodel.py`, which builds Cantera reactions in memory. Own ticket.
- **`rmgpy/rmg/output.py:647`** — sets `duplicate = True` for every `Multi*` in
  `save_diff_html`. Same false premise as HIGH 3, but it is the `diffModels.py` HTML display
  path and emits no mechanism. Enumerated, not repaired.
- **`test_family_still_generates_the_sei_reaction`** (`preflightDeckFamilyExclusionTest.py`)
  — **status changed, and downgraded.** Round 74 reported it hanging past 900 s (two runs,
  `timeout 900` → exit 124). In this round's full top-level run it **PASSED**, inside a
  suite that completed in 2097 s. So the hang is not deterministic, and "it hangs" is not a
  safe claim. What can be said: it is slow enough to have exceeded a 900 s budget twice, and
  it is `@pytest.mark.database` so `make test` excludes it. It has still never been run
  against the base, so "pre-existing" remains an inference and not a measurement.

## 10. One incidental observation, not repaired

`Lindemann` and `Troe` are **not** subclasses of `ThirdBody` in RMG (all three derive
directly from `PDepKineticsModel`). So the `and not isinstance(kinetics, (Lindemann, Troe))`
guard inside `chemkin_third_body_token_shape`, added in round 74, is dead. It is harmless and
protective if the hierarchy ever changes, so it stays — but it reads as load-bearing and is
not, and that is worth knowing before anyone reasons from it.
