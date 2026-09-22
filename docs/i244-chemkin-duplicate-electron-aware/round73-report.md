# Round 73 — the key, the scope, and the generator

Base engine for every before-arm: `/home/alon/Code/RMG-Py-plasma` at `311818121`, built.
Branch engine: this worktree, `4172ed954` before the repair.
Interpreter: `/home/alon/anaconda3/envs/rmg_env/bin/python`. Cantera: `/home/alon/Code/Cantera/build/python`.

Every log named below is under `docs/i244-chemkin-duplicate-electron-aware/logs/`.

## 1. The three regressions, reproduced here on both arms before anything changed

Script: `round73_repro.py`, written from scratch for this round; self-contained, so it runs
against any checkout. It names the loaded engine by a value, not a path.

| | BASE `311818121` | BRANCH `4172ed954` | AFTER |
|---|---|---|---|
| **D1** permuted pre-marked pair | 2 `DUPLICATE`, **`cantera.Solution` LOADED** | 0 `DUPLICATE`, **`cantera.Solution` REJECTED** | 2 `DUPLICATE`, **LOADED** |
| **D2** core save then core+edge save | flags `(False, False)`; Cantera entry `duplicate=False` | flags `(True, True)`; **lone `duplicate: true` shipped** | `(False, False)`; `duplicate=False` |
| **D3** generator argument | `TypeError` — refuses loudly | 1 reaction from a list, **0 from a generator** | 1 and 1 |

Logs: `R73-repro-BASE-*`, `R73-repro-BRANCH-*`, `R73-repro-FIXED-*`.

**The near-miss the brief warned about is real and is in the log.** On the branch arm D1 reads
`ck2yaml.convert_mech ACCEPTED; cantera.Solution REJECTED`. Conversion is not validation. The
repro and every new test load the mechanism.

Two notes on the evidence, so it reads straight:

- `R73-repro-BRANCH-stdout.log` was produced by the first revision of the script. Its banner
  probed `mark_duplicate_reactions.__doc__` for "group-level", which is how the brief said to
  identify the engine — correct at `4172ed954`, and wrong after the repair, which moves the
  group-level authority to `chemkin_duplicate_flags`. The banner now asks both functions and
  names whichever answers. The BASE and FIXED logs are from the current revision; the BRANCH
  log is not, and nothing but its banner line differs.
- D1's "flags after the writer" line reads `(True, True)` on the repaired arm and is not
  evidence of marking: the repaired renderer does not touch the flags at all, so what the
  line shows is the input surviving. The `DUPLICATE` count and the Cantera verdict are the
  measurement.

## 2. Defect 1 — the key is now a multiset per side

`chemkin_duplicate_group_key` sorts each side's `id()` tuple. Everything else in the key is
unchanged: class, `specific_collider`, per-side electron placement, pressure dependence,
reversibility, and the reversible case's unordered pair of sides.

**Why a multiset and not the rendered equation.** The rendered string is what Cantera reads, so
keying on it would also catch two distinct `Species` objects that render to the same identifier.
It costs three things: `mark_duplicate_reactions` is not given the deck's species list, so the
key cannot build the identifier; `get_species_identifier` can raise, and a duplicate check is a
bad place to raise from; and folding two `Species` the rest of RMG holds apart would hide a
deck-level defect rather than fix it. The multiset gets the whole of the reported defect for
none of that.

The docstring's old defence — that permutations "were not duplicates before and are not now;
that is a separate question" — is retracted in the docstring itself, with the reason: it stopped
being a separate question when the function gained the power to *clear*.

### The census, with its denominator

`round73_library_census.py`, log `R73-library-census-stdout.log`.

```
libraries loaded              : 190 of 190
reactions examined            : 35277
entries whose answer CHANGES  : 4
groups merged by the multiset : 2
label/formula collisions      : 0
```

The unit is one library, because that is the unit in which a permuted pair can be authored.
Species identity inside a library is taken from the species label, since freshly loaded entries
do not share `Species` objects; zero labels resolved to more than one molecular formula, so the
substitution held everywhere.

The four are two pairs, in `NIST_Fluorine/full` and `NIST_Fluorine/seed`:

```
CHF2-CF3 + C3H7 <=> CF3-CF2 + C3H8
C3H8 + CF3-CF2 <=> C3H7 + CHF2-CF3
```

One reversible entry written twice, sides swapped and each side permuted. Checked directly
(`R73-census-fluorine-detail-stdout.log`): **both arrive `duplicate=True`** — RMG's own
`KineticsLibrary.check_for_duplicates` marks them at load. Under the ordered key each was a
singleton and both flags were cleared, so any deck containing that library shipped one equation
twice with no `DUPLICATE` line. Neutral fluorine chemistry; no electron anywhere in it.

## 3. Defect 2 — the flag stops being the channel

Of the two repairs the brief offered, this is the second: not "make the recompute not outlive
the render" but "stop the flag being how the render is told".

- `chemkin_duplicate_flags(reactions)` computes the group answer and returns a list of
  booleans. It mutates nothing.
- `write_kinetics_entry` takes `duplicate=None`, falling back to `reaction.duplicate` for
  callers holding one reaction and no list.
- `render_chemkin_file` and `render_chemkin_surface_file` compute the answer for their own list
  and pass it per entry.
- `save_chemkin` no longer recomputes centrally, and no longer passes
  `check_for_duplicates=False` to its renders.
- `mark_duplicate_reactions` survives as the mutating wrapper — long-standing public API, and
  `rmgpy.rmg.model.mark_chemkin_duplicates` genuinely wants the flags left on the objects.

**What it costs.** Each file keys itself, so a verbose render repeats the work of its plain
render: two linear passes per file instead of one for the whole save. The alternative was an
extra parameter threading a precomputed answer from `save_chemkin` into each render, which is
the same channel problem one level up. It also *fixes* a second scoping error nobody reported:
`save_chemkin` used to key the gas file's entries against the surface file's.

**Every reader of `Reaction.duplicate`, named.** Writers first:

| site | reads the flag? | consults the group answer? |
|---|---|---|
| `rmgpy/chemkin.pyx` `write_kinetics_entry` | only when `duplicate=None` | yes, via its renderer |
| `rmgpy/chemkin.pyx` `render_chemkin_file` / `render_chemkin_surface_file` | when `check_for_duplicates=False` | yes |
| `rmgpy/yaml_cantera2.py:681` | **yes, directly** | **no** |
| `rmgpy/reaction.py:570` `to_cantera` | **yes, directly** | **no** |
| `arkane/pdep.py:506` | via `write_kinetics_entry(duplicate=None)` | **no** |
| `rmgpy/chemkin.pyx:301` | yes — but this is the *reader* parsing a deck back in, not a writer |

`rmgpy/rmg/main.py` is a reader only indirectly: `generate_end_of_run_cantera_files` runs the
Cantera writers on the core after `save_chemkin_files` has written core and core+edge. That
sequence is the whole of Defect 2.

**The Arkane claim, checked.** Half right.

- `arkane/pdep.py` — **confirmed**. It imports `write_kinetics_entry` and calls it per network
  reaction (line 506), never passing through any renderer. It cannot consult a group answer
  because it has no list to key over; its local `duplicate` variable is an unrelated
  "comment this one out" flag for a repeated reverse channel.
- `arkane/kinetics.py` — **refuted**. Its `write_chemkin` (line 290) formats its own line with
  `'{0!s:51} {1:9.3e} ...'`, appends it to `chem.inp`, and never reads `reaction.duplicate` or
  emits a `DUPLICATE` line at all. It is not a duplicate-flag writer.

**Not repaired, reported.** `yaml_cantera2` and `Reaction.to_cantera` read the flag with no
group recompute of their own, so a Cantera artifact still carries whatever
`rmgpy.rmg.model`'s incremental pairwise sweep left. That is base behaviour, predates this
branch, and is outside this ticket; what this repair guarantees is only that the Chemkin
writers no longer *feed* them a wrong answer. Worth its own ticket.

## 4. Defect 3 — materialized, and nothing upstream wanted laziness

`render_chemkin_file` and `render_chemkin_surface_file` do `reactions = list(reactions)` on
entry. Checked before materializing — every caller in the tree already passes a sequence:

- `rmgpy/chemkin.pyx` `save_chemkin`: `core.reactions + edge.reactions`, a new list;
- `rmgpy/tools/mergemodels.py:114`, `scripts/standardizeModelSpeciesNames.py:110`,
  `scripts/thermoEstimator.py:85`: lists;
- every test: lists.

And the base engine already required a sized sequence — its `mark_duplicate_reactions` called
`len(reactions)`, which is the `TypeError` the base arm raises. So materializing costs one
pointer array that the caller was already holding, and it restores the contract
`mark_duplicate_reaction` has documented for years.

## 5. Defect 4 — the repointed test, and the claim checked

`test_a_cross_group_comparison_does_not_clear_another_groups_flags` had all four inputs arriving
`True` and all four expectations `True`. `undeclared_b` now arrives `False` and must be *set*.
The defect under test is untouched: the cross-group contamination runs between `declared_a`,
`declared_b` and `undeclared_a`, all three still pre-marked.

Measured against a no-op authority (`chemkin_duplicate_flags` returns the flags it was handed),
log `RED-D4-repoint-comparison-stdout.log`:

```
--- the version at HEAD (4172ed954), all four inputs pre-marked True ---
1 passed
--- the repointed version, one input arriving False ---
1 failed
```

The brief's claim holds exactly.

## 6. Every new and repaired check, shown red

`round73_red_arms.py` breaks one thing per arm, rebuilds, runs the file, restores, rebuilds.
`round73_arm_discriminators.py` prints four values off the *built extension* so that a build
that silently did not take shows up as the wrong line moving.

| arm | discriminator that flips | tests RED |
|---|---|---|
| `D1-ordered-key` | D1 permuted pair keys alike → False | `test_a_permuted_pair_is_not_cleared_into_a_deck_cantera_rejects` |
| `D2-leaking-recompute` | D2 render leaves flags untouched → False | `test_a_render_leaves_every_reactions_flag_exactly_as_it_found_it`, `test_a_core_plus_edge_save_does_not_mark_a_core_only_cantera_entry` |
| `D3-consumed-generator` | D3 generator survives the render → False | `test_a_generator_of_reactions_is_written_not_consumed` |
| `D4-noop-authority` | D4 authority overrides its input → False | 14, including the repointed cross-group test |

In each of D1–D3 exactly one discriminator moved and it was the arm's own. Restored engine:
all four `True`, 24 passed (`GREEN-discriminators-stdout.log`, `GREEN-i244-suite-stdout.log`).

**A false start, recorded because it nearly shipped.** The first revision of the harness
restored the source with `git checkout -- rmgpy/chemkin.pyx`. That restores from the *index*,
and the repairs were not staged — so it reverted them to `4172ed954`. Arms D2, D3 and D4 then
all measured the unrepaired engine and produced the same four failures while appearing to
measure three different breakages. The harness now keeps a byte backup and restores from that,
with `shutil.copy` plus `os.utime` rather than `copy2`, because `copy2` preserves the backup's
mtime and Cython then skips the rebuild — which is how the *first* corrected run still left the
broken D4 engine loaded after its restore. Both traps are commented in the harness. The
discriminator file exists because of the first one.

## 7. Verifier

1. Three regressions reproduced here on both arms before any change, D1 through
   `cantera.Solution` — §1. ✅
2. `examples/rmg/minimal`, run end to end on each engine, edge decks through
   `cantera.Solution` (`R73-minimal-cantera-validation-stdout.log`):

   ```
   min-BASE /chemkin/chem.inp        DUPLICATE=0   LOADED  (26 species, 66 reactions)
   min-BASE /chemkin/chem_edge.inp   DUPLICATE=11  REJECTED: Kinetics::checkDuplicates
   min-FIXED/chemkin/chem.inp        DUPLICATE=0   LOADED  (26 species, 66 reactions)
   min-FIXED/chemkin/chem_edge.inp   DUPLICATE=10  LOADED  (172 species, 398 reactions)
   ```

   Exactly one `DUPLICATE` line fewer, 398 reactions either way. The headline stands and the
   ordered-pair deck is valid too (§1, D1). ✅
3. Key basis chosen and argued, census with denominator — §2. ✅
4. Scoping repair argued, every reader named, Arkane claim checked both ways — §3. ✅
5. Generator fixed, laziness checked — §4. ✅
6. The repointed test shown red against a no-op authority — §5, §6. ✅
7. Docstrings corrected: the permutation claim, "what every writer calls", "one pass". ✅
8. `git diff --check` clean; `git status` clean apart from `docs/pm-verify/` and the
   `/dev/null` dotfiles. ✅

## 8. Suite

`chemkinTest`, `i244ChemkinDuplicateElectronTest`, `i134DuplicateElectronsTest`,
`i126ChemkinElectronOrderTest`, `i167CanteraExportPathTest`, `i135TdepRoundTripTest`,
`yaml_cantera1Test`, `yaml_cantera2Test`, `reactionTest`, `plasmaExportTest`,
`rmg/modelTest` — `R73-suite-chemkin-export-stdout.log`:

```
5 failed, 564 passed, 2 skipped, 1 xfailed
```

The five are the pre-existing `i134` failures of §9; their node list is byte-identical to the
base engine's, diffed.

## 9. Findings outside the repair

**`test/rmgpy/i134DuplicateElectronsTest.py` has 5 failures that are not this branch's.**
Identical node list on `311818121` and here, in isolation on both
(`R73-i134-BASE-stdout.log` against `R73-suite-chemkin-export-stdout.log`). The visible cause is
database content: a library the test expects to hold one lithium entry now also holds
`[Arp] => [Ar]`, so `assert len(reactions) == 1` sees 2. Pre-existing, database-side, and left
alone here.

**`.gitattributes` needed a second lever.** `-whitespace` does not suppress `git diff --check`'s
leftover-conflict-marker report; that is keyed on the marker *length*. Verified in a scratch
repo: with `-whitespace` alone a log containing `<<<<<<< HEAD` is still reported, with
`conflict-marker-size=200` it is not. Both committed logs now pass without being edited.

**"One pass, linear" was wrong and is now stated true**: two passes, and linear in reactions
plus participant occurrences, since each key sorts its own two sides.
