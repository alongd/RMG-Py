# I-244 round 74 — the key stops carrying class, and the Cantera writer keys its own list

Rework of the two HIGHs and two LOWs returned by the round-73 adversarial review. Base of
this round: `242431e00`. Every deck verdict below comes from `cantera.Solution`, never from
deck text and never from `ck2yaml.convert_mech` alone.

## 1. HIGH 1 — cross-class identity

**Defect.** `chemkin_duplicate_group_key` carried `reaction.__class__`. Class is absent from
Chemkin serialization, so a `LibraryReaction` and a `TemplateReaction` over the same
participants render one identical equation — but keyed apart, they became two singleton
groups and the singleton path *cleared* both. The deck went out with zero `DUPLICATE` lines
and `Kinetics::checkDuplicates` rejected it. While doing it, the module logged

```
Marked reaction H(1) + OH(2) => H2(3) + O(4) as not duplicate because no other reaction writes its Chemkin equation.
Marked reaction H(1) + OH(2) => H2(3) + O(4) as not duplicate because no other reaction writes its Chemkin equation.
```

— the same sentence twice, for the same equation, each time asserting the very thing the key
had failed to notice. Captured at `logs/R74-RED-both-highs-stdout.log:134-135`.

**Provenance, which changes the argument.** The term is not an invention of this branch. It
was carried verbatim from stock's pairwise predicate, `git show 311818121:rmgpy/chemkin.pyx`:

```python
if reaction1.__class__ != reaction2.__class__:
    # TemplateReaction, LibraryReaction, and PDepReaction cannot be duplicates of one another.
    # RHW question: why can't TemplateReaction be duplicate of LibraryReaction, in Chemkin terms?
    #               I guess it shouldn't happen in RMG.
    continue
```

Stock's own author questioned the rule inline and never answered it. More decisively, **the
term does not survive the move from a predicate to a group key**: in the predicate a class
mismatch means `continue`, so the pair is never marked *and never cleared*; in a group key
the same term splits the pair into singletons the clearing path then strips. One term,
opposite consequences. This is the round-73 lesson recurring — clearing is the dangerous
direction, and every argument that justified a predicate's narrowness inverts when the
function gains the power to clear.

**Repair.** The class term is removed from the key, and the docstring records the retraction,
the stock provenance, and the direction argument.

**Was any cross-class regression justifying it?** No, checked three ways rather than assumed.
No test in the tree asserts cross-class non-duplication. Removing a term from a key can only
*merge* groups, never split them, so it can only add `DUPLICATE` lines and never remove one;
every property that changes the rendered equation is still keyed, so the entries it merges are
exactly the entries that render alike. And measured on a real mechanism: `examples/rmg/minimal`
re-run with the class term gone produces `chem.inp`, `chem_edge.inp` and
`species_dictionary.txt` **byte-identical** to the round-73 run. Zero change on ethane
pyrolysis; the term was inert there and lethal on the constructed case.

## 2. HIGH 2 — the scope leak survives through production marking

**Defect.** `rmgpy.rmg.model` marks duplicates on the shared core and edge reaction objects
while the model grows — `model.py:838`, `:842`, `:1967` — before any writer runs. The repaired
Chemkin renderers key their own list and deliberately leave the objects' flags alone, which is
exactly what lets those production marks survive the render. `rmgpy.yaml_cantera2` then read
`Reaction.duplicate` directly, so a core-only export inherited the core+edge answer: a lone
`duplicate: true`, rejected by `Kinetics::checkDuplicates`.

Reproduced at `logs/R74-RED-both-highs-stdout.log`:
`No duplicate found for declared duplicate reaction number 0 (H(1) + OH(2) => H2(3) + O(4))`.

**Repair.** `_collect_reaction_entries` materializes its list, computes
`chemkin_duplicate_flags` over it, and passes the answer per entry;
`reaction_to_dict_list` takes `duplicate=None` and falls back to the reaction's own flag only
for single-reaction callers. Gas and surface entries are collected by separate calls, so each
keys only against its own phase. The `MultiArrhenius`/`MultiPDepArrhenius` expansion passes
`duplicate=True` explicitly, since those sub-entries all carry one equation by construction —
the same rule the Chemkin MULTI block follows.

**Production reach, measured.** `CanteraWriter2` is a real writer, attached when
`generateCanteraYAML2=True`, and `main.py:1495` calls `cantera2/chem.yaml` "the authoritative
Cantera artifact for this run" and the only writer that can represent a plasma mechanism. It
calls `save_cantera_model(rmg.reaction_model.core, ...)` — the core alone — at every
iteration, on the objects model growth has marked. The structural condition for the leak is
therefore met on every such run.

**But it does not fire on `examples/rmg/minimal`, and that is worth stating plainly.** Both
arms were run end-to-end with `generateCanteraYAML2=True`: 23 core exports and 26 edge exports
each. With the repair, all 23 core exports load. With the leak deliberately reinstated, all 23
core exports still load, and every one of them is **byte-identical** to the repaired run's.
The reason is that on ethane pyrolysis no core reaction is ever marked duplicate at all — the
ten marks in the edge deck are edge-edge pairs — so there is no core+edge answer to leak onto
a core entry. The defect is real and the repair is right, but its blast radius on this
mechanism is zero, and a reader should not be told otherwise. It needs a mechanism in which a
core reaction's only mate lives on the edge.

## 3. The two tests the review said could not fail

- `test_a_render_leaves_every_reactions_flag_exactly_as_it_found_it` asserted only that
  nothing changed, which a renderer that did nothing would satisfy — all three flags arrive
  `True` and the expectation was that all three stay `True`. It now pins the rendered deck
  first (`[True, True, False]`, the lone entry cleared) and the flags second, so the answer
  written and the flags left alone are required to *differ*. Confirmed: under the no-op
  authority arm this test now fails, where before it passed.
- `test_a_core_plus_edge_save_does_not_mark_a_core_only_cantera_entry` built both reactions
  unmarked and closed by calling `reaction_to_dict_list` on objects nothing had ever marked,
  which no implementation could fail. It is renamed
  `test_a_core_plus_edge_save_does_not_store_its_answer_on_the_reactions`, its docstring now
  states its scope and names what it cannot prove, and the Cantera claim moves to the new
  class that reaches the production path.

New class `TestProductionMarkingDoesNotReachTheCanteraWriter` starts from production marking,
not from fresh objects, and ends at `cantera.Solution`. It carries a tie-back asserting
`rmgpy.rmg.model` still contains all three marking calls, so it fails loudly rather than
vacuously if model growth stops marking.

## 3b. A third defect, found by fixing HIGH 2 — and it was mine

Applying the group key to the Cantera writer's list pointed it at reactions the Chemkin
renderers had never handed it, and `plasmaExportTest.py` went red. The cause was a real defect
in the key, introduced in round 73 and invisible until now.

`chemkin_duplicate_group_key` carried pressure dependence as a **boolean**. The rendered
equation does not: `chemkin.pyx:1976-1983` writes `+M` for a `ThirdBody`, `(+M)` (or
`(+collider)`) for a `Lindemann`/`Troe`, and nothing at all for a `PDepArrhenius`. A boolean is
true for all three, so a `ThirdBody` writing `A+B+M=>C+M` was grouped with a `Troe` writing
`A+B(+M)=>C(+M)` — two different equations in one group, which with the group able to set flags
marks the three-body as a duplicate of an equation nothing else writes. Chemkin and Cantera
refuse a lone `DUPLICATE` exactly as they refuse an undeclared pair.

Repaired by keying on the *shape* of the third-body token. The shape now comes from one
function, `chemkin_third_body_token_shape`, which `write_kinetics_entry` also uses to build the
token — so the key and the writer cannot drift apart about what reaches the equation, which is
how this defect existed at all. The key keeps the pressure-dependence boolean alongside it, so
the long-standing mixed-pressure-dependence behaviour that round 73's tests pin is unchanged.

**And the fixture it broke was itself invalid.** `test_gas_export_is_byte_identical_to_the_base_commit`
compares against a committed golden, `test/rmgpy/test_data/plasma_export/t0_one_gas.yaml`. That
golden holds three entries writing `O(2) + O2(3) => O3(4)` with only two of them marked.
`cantera.Solution` **rejects** it: *"Undeclared duplicate reactions detected"*. A byte-identity
test could never see that — it compares text against text, and both sides were equally invalid.

The goldens are regenerated. The only content difference is added `duplicate: true` lines,
three in the gas golden and two in the surface one, plus the `generator:` line the test already
normalizes; no rate constant, reference temperature or thermo coefficient moved, so the
`T0 = 1 K` invariant those tests exist to guard is untouched. Both regenerated goldens load.
A new test, `test_the_goldens_are_mechanisms_cantera_can_load`, is the tripwire: it fails
against the golden as committed at `242431e00` and passes against the regenerated one, so a
future regeneration cannot quietly reinstate an invalid deck.

This is the third time in two rounds that the answer was wrong in a way only *loading* the
mechanism could reveal — first stock's edge deck, then the cross-class clear, now a test
fixture. Text comparison has now missed it three times for three different reasons.

## 4. RED/GREEN arms

`round74_red_arms.py`, restoring from a byte backup with `shutil.copy` + `os.utime` — never
`git checkout`, never `copy2`; both traps are commented in the source and both have cost this
campaign a round. Each arm prints four discriminators read off the built extension.

| arm | discriminator moved | tests failed |
|---|---|---|
| `D1-class-back-in-key` | `key_ignores_class` → false | 1 — the cross-class test only |
| `D2-cantera-reads-the-flag` | `cantera_writer_keys_its_list` → false | 1 — the production test only |
| `D3-third-body-bool-only` | `key_separates_third_body_from_falloff` → false | 1 — the three-body test only |
| `D4-noop-authority` | `flags_recompute_clears_a_lone_mark` → false | 19, incl. all three new tests and the repaired nonmutation test |
| restored | all five true | 0 — 29 passed |

Each arm moved exactly one discriminator and failed exactly its own guard, which is what
distinguishes four real measurements from four runs of the same broken build.

Separately, the golden tripwire was shown red by loading the golden as committed at
`242431e00` (`REJECTED: Kinetics::checkDuplicates — Undeclared duplicate reactions detected`)
and green against the regenerated one.

## 5. The two LOWs

**The complexity claim** at `chemkin_duplicate_flags` said "linear in the number of reactions
plus the number of participant occurrences, since each key sorts its own two sides" — a
sentence whose subordinate clause contradicts its main claim. Corrected: the cost is
`sum over reactions of k log k` in the participants per side. It *is* linear in the number of
reactions, which is the term the rewrite was about.

**`.gitattributes`** — `conflict-marker-size=200` is reverted, not argued. Reverting it
surfaced something the suppression had hidden: the check reports **three** logs, not the two
known ones. `R73-suite-chemkin-export-stdout.log:752` accumulated while the detector was off
and nobody saw it — precisely the cost of silencing a detector repo-wide to quiet two known
lines.

All three are false positives, and are recorded rather than suppressed. The flagged lines are
pytest's short-summary separator, which at this terminal width begins with exactly seven `=`
followed by a space — git's pattern for a conflict's middle marker. (The reverted comment
described them as seven `<`; that was wrong, and is corrected.) `-whitespace` is kept, since
the logs are evidence and must stay byte-identical; the honest state is that
`git diff --check 311818121..HEAD` reports three lines on this branch.

## 6. No-regression

| deck | DUPLICATE | `cantera.Solution` |
|---|---|---|
| minimal core, round 74 | 0 | LOADED — 26 species, 66 reactions |
| minimal edge, round 74 | 10 | LOADED — 172 species, 398 reactions |

Byte-identical to the round-73 decks — `chem.inp`, `chem_edge.inp` and
`species_dictionary.txt` all three, re-run after the final engine including the third-body key
refinement. The measurement this branch exists for is unmoved: base `311818121` still emits an
edge deck `Kinetics::checkDuplicates` rejects (`DUPLICATE=11`), and this branch still emits one
Cantera accepts (`DUPLICATE=10`, 172 species / 398 reactions).

Suites, all with the engine as committed: `test/rmgpy/*.py` 995 passed; `test/rmgpy/rmg/` plus
`test/arkane/` 356 passed; `test/rmgpy/data/` 358 passed. Zero failures anywhere. The i134 five
that were red in round 73 are gone — I-255 landed the repair, and this branch did not touch
that file.

## 7. One thing left open, reported rather than fixed

`test/rmgpy/preflightDeckFamilyExclusionTest.py::ElectrochemicalFamilyPreservedTest::test_family_still_generates_the_sei_reaction`
**hangs** — it ran past a 900 s timeout twice. It is `@pytest.mark.database`, so `make test`
excludes it and round 73's suite deselected it along with 120 others; it was never green on this
branch and has never been run by this campaign's suites. It loads a kinetics family and calls
`generate_reactions_from_families`, touching no Chemkin writer and no Cantera writer, so nothing
changed in rounds 73 or 74 is on its path. **I did not verify it against the base engine** — that
needs a base rebuild — so "pre-existing" is an inference from the call path, not a measurement,
and it is recorded here as such rather than claimed.
