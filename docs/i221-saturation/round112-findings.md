# Round 112 — the map is keyed on the species, and the reaction is made of occurrences

Base `i221-saturation-no-atomtype@71ae97bd5`. Commits `4f74c023f` (HIGH 1, 2), `e3f646b2d`
(HIGH 3), `c44390545` (HIGH 4, 5), `898b85d8d` (MEDIUM). All local; nothing pushed, merged or
rebased.

All five HIGHs and the MEDIUM are reproduced red and closed. The end-to-end gate passes for the
half this branch can reach; the other half ("a build carrying the other campaign work") needs a
merge this round may not make, and is reported as **not reached** (§7).

---

## 1. HIGH 1 and HIGH 2 — pairs are occurrences, resolved in one place

The fix is in the representation, not at the call sites. `reaction.pair_occurrences(pairs,
reactants, products)` turns pairs into side-specific `(reactant position, product position)`
tuples: the k-th time a species is named on a side means its k-th occurrence there, and a member
that is not on its side raises `ValueError`. `pairs_from_occurrences` rebuilds pairs on new lists
from those positions. Every remapper goes through the pair: `Reaction.copy`, `model.py`'s
`make_new_reaction` (which now carries positions instead of calling `forward.reactants.index(...)`
/ `products.index(...)`), and the merge/ensure paths.

`test/rmgpy/i221OccurrencePairsTest.py`, 18 tests, asserting on positions and object identity, not
on labels or species sets. At `71ae97bd5` (`logs/round112-tests-red-at-71ae97bd5.log`):
**13 failed, 5 passed**. The five that pass at base are the controls, and they must pass there:
species on both sides, a thermal reaction, no pairs, reverse pairs, and a thermal model build. The
brief's five required controls are all present:

| control | test outcome |
|---|---|
| species twice on one side | positions preserved across copy/model |
| species on both sides | preserved (passes at base too) |
| duplicated electron products | `(0,1)`/`(0,2)` stay distinct, not collapsed to the last copy |
| missing pair member | **raises** `ValueError` |
| ordinary thermal reaction | unchanged (passes at base too) |

Focused suite at the HIGH 1/2 commit (`logs/round112-focused-high12.log`): **672 passed, 4
skipped, 52 deselected**.

### The `.index` audit (HIGH 2)

No `list.index()` call is left on a pair path in the gated files. Every one that remains searches
a list of distinct objects, so it cannot match the wrong occurrence:

| site | list searched | safe because |
|---|---|---|
| `rmgpy/data/kinetics/isotopes.py:310` | `mol.atoms` | atoms are distinct objects |
| `family.py:3923/3924/4024/4025` | atom indices built for the tree | distinct atoms |
| `family.py:3949` | `vals.index(min(vals))` | numbers, not species |
| `family.py:4340` | an integer inverse permutation | a permutation has no repeats |
| `rmgpy/rmg/pdep.py:409/413/444/450/469` | `isomer_spcs`, keyed on unimolecular `reactants[0]`/`products[0]`/`source[0]` | isomers are unique and each key is a single species |
| `pdep.py:634` | `networks` | distinct networks |
| `pdep.py:934` | `configurations` | distinct configurations |
| `model.py:1526` | `prunable_species` | a deduplicated species list |

The four other hits are docstrings or comments: `reaction.py:2342`, `library.py:165`,
`family.py:693` and `model.py:653`.

## 2. HIGH 3 — a transport cannot exist without the complete reducers

This is the mechanical option. `test/rmgpy/i221TransportCensusTest.py` walks the source tree for
every process transport: `Pool`, `Process`, `Pipe`, `Queue`, `concurrent.futures`, and pickler
construction. It fails for any module that holds one but does not call `install_complete_reducers`
at import. `rmgpy/qm/main.py` now installs them at module level.

At the parent of the fix (`logs/round112-transport-census-red.log`): **3 failed, 15 passed**. At
the fix (`…-green.log`): **19 passed**.

**Census, extended rather than re-derived:**
- Process transports: **5 calls in 3 modules**. They are `rmg/react.py:80` `Pool`,
  `qm/main.py:266` `Pool`, and in `family.py` the `Pool` at 4343 plus the `Pipe`/`Process` pair
  at 5495/5497. The `psutil.Process` calls and `pathfinder.py`'s FIFO `Queue` are not transports.
- Modules that import a pickler: **1**, `family.py`.
- Pickle writes to storage: **0**.
- `__reduce__` exposure that the reducers do not reach: `Reaction:322`, `LibraryReaction:119`,
  `ReactionModel:162`. This is accepted because no transport outside `family.py` pickles them, and
  `family.py` goes through `_CompletePickler`.

## 3. HIGH 4 and HIGH 5 — enumeration is fresh, and "could not look" is not "found nothing"

`_enumerate_quarantines` now returns `(found, incomplete)`. A loaded family is re-read through
`resolve_quarantine`, the signature-validated disk path, instead of being served from its in-memory
copy, so a quarantine rewritten on disk after load is seen (HIGH 4). If the enumeration could not
look at a place, it records the gap. The caller yields everything it did find and then raises
`QuarantineEnumerationError`, so an unreadable families directory can no longer read as zero
quarantines (HIGH 5).

At the parent (`logs/round112-quarantine-red.log`): **6 of 7 new tests fail**. At the fix
(`…-green.log`): **361 passed, 2 skipped, 36 deselected**. Two existing tests that pinned the
silent behaviour were updated.

**Catch-and-return-empty census:** 6 sites closed.
1. enumeration `listdir`
2. enumeration `lexists`
3. the silent skip of an unreadable loaded manifest
4. has-any `stat`
5. has-any `listdir`
6. has-any `lexists`

The remaining `except` sites are left because each one is either distinguishable from "nothing
found" or cannot affect a verdict:
- the `UNREADABLE` sentinel;
- the path where realpath returns `None`, which is reported as unanswered;
- `_manifest_signature` returning `'unreadable'`;
- the AST parse returning an explicit `False`;
- a cosmetic `close` and a `repr`;
- the wording in `_why_unanswered`.

## 4. MEDIUM — resolved in the refusing direction

`_CompletePickler.reducer_override` refuses, at pickling time and by name, any subclass of a class
that has a lossy reducer when the subclass itself is not registered. One subclass is accepted on
purpose and the reason is recorded: `_ACCEPTED_UNREGISTERED_SUBCLASSES = {'arkane.encorr.data.Molecule': …}`.

The private Cython memos (`_fingerprint`, `_inchi`, `_smiles` on `Molecule` and `Species`) are
documented as accepted. They are caches that get recomputed, and a test pins them against the
`.pxd`, so adding a new private attribute fails loudly.

The old whole-pickle test only compared scalar fields. The new one compares nested state, so a loss
inside a nested object now fails it.

One case cannot be caught: an unrelated class that is itself lossy through `__reduce__`, with no
registered base class. There is nothing to detect it from, so it is accepted.

At the parent (`logs/round112-medium-red.log`): **4 failed, 4 passed**. At the fix
(`…-medium-green.log`): **519 passed, 4 skipped, 36 deselected**.

## 5. End-to-end gate (verifier 10, first half) — PASS

Probe: `docs/i221-saturation/probes/round112_end_to_end_probe.py`. Log: `logs/round112-e2e.log`.
Result: **49 ok, 0 failures, exit 0**.

It runs the real `Plasma_Electron_Impact_Ionization` family from `RMG-database-plasma@96f2afa4a` on
three reactants: Li, CH3, and metastable Ar (`1 Ar u2 p3 c0`). The stages are:
1. generate
2. `resolve_electron_placement`, then `generate_pairs`
3. `copy()`, checking that `pair_occurrences` is the same before and after and that nothing is aliased
4. `saturate_for_estimation` on every radical
5. `update_atomtypes`
6. canonical `copy()`
7. `CoreEdgeReactionModel.make_new_reaction`, checking that the pair positions equal the canonical
   ones and that the pair members are the model's own species
8. thermo

A bare `AtomTypeError` at any stage counts as a failure. There were none.
- Metastable Ar saturation and thermo are refused by name with `SaturatedStructureError`.
- Ar⁺ thermo raises a loud `DatabaseError` (no group node). The argon cation library was
  deliberately not loaded; production loads it.
- The view's pair positions stay distinct through copy: `[(1,1),(0,0),(0,2)]` for Li and Ar, and
  `[(0,0),(1,1),(1,2)]` for CH3.

## 6. Suites

| suite | result | log |
|---|---|---|
| RMG-Py, whole tree minus `data/rmgTest.py` | 5 failed, 3820 passed, 41 skipped, 1 xfailed (1:28:59) | `round112-full-suite.log` |
| RMG-Py, `test/rmgpy/data/rmgTest.py` alone | 2 passed | `round112-full-suite-rmgTest.log` |
| `test/rmgpy/i134DuplicateElectronsTest.py`, tip, isolated | 5 failed, 178 passed, 1 xfailed | `round112-i134-tip.log` |
| the same at `71ae97bd5` | 5 failed, 178 passed, 1 xfailed, **the same five** | `round112-i134-at-71ae97bd5.log` |
| RMG-database-plasma `test/` at the tip | 1 failed, 248 passed | `round112-database-suite.log` |
| the same at `71ae97bd5` | 1 failed, 248 passed, **the same test** | `round112-database-suite-at-71ae97bd5.log` |

The whole tree cannot be collected in one run because two test files share a basename:
`test/rmgpy/data/rmgTest.py` and `test/rmgpy/rmg/rmgTest.py`. `--import-mode=importlib` avoids that
but breaks the yaml_cantera tests, which import their sibling `cantera_yaml_comparer`. So the tree
runs in two halves.

**All 5 RMG-Py failures are in `i134DuplicateElectronsTest`, and they fail the same way at
base.** Each one asserts `PlasmaRadiativeRecombination has 2 entries, expected 1`. That count is
the database's, not this branch's. The fixture pins a library that grew when the database checkout
moved, so the failure is database drift. It needs a re-pin on the I-134 side and is outside these
gates. The base run needs an `rmgrc` with an absolute path. With the template's relative path, the
worktree at `/tmp` resolves a nonexistent database and gives a different, meaningless 6F+5E.

**The database failure is the same at base.** It is caused by I-221's core change (`73c8214f7`),
not by this round.
`test_argon_cation_buildtime.py::test_the_dimer_cation_still_builds_but_is_now_refused_one_layer_earlier`
pins the bare `AtomTypeError`. I-221 now raises `SaturatedStructureError` from it, deliberately, with
the `AtomTypeError` chained as the cause. The test needs to be re-pinned to `SaturatedStructureError`
when I-221 merges. That is a cross-branch edit to the database repo, outside this round's gates.

Round-110 database counts are not comparable, because the database checkout has since moved to
`96f2afa4a`.

## 7. Not reached

Verifier 10, second half: *"a build carrying the other campaign work that no longer dies at
`saturate_radicals` with the bare `AtomTypeError`"*. Producing that build means merging this branch
with the rest of the campaign, and this round may not merge. The database failure in §6 is the one
concrete piece of evidence available: the pinned database test from the other side of that merge
now receives `SaturatedStructureError` instead of the bare `AtomTypeError`, which is the behaviour
the gate asks for.

## 8. Commands

All runs use `PATH=/home/alon/anaconda3/envs/rmg_env/bin:$PATH` and capture both streams with
`> >(tee -a $L/<name>.log) 2> >(tee -a $L/<name>.stderr.log >&2)`, where
`L=docs/i221-saturation/logs`.

```bash
# focused (MEDIUM green)
python -m pytest test/rmgpy/data/kinetics/ test/rmgpy/i221OccurrencePairsTest.py \
  test/rmgpy/i221TransportCensusTest.py test/rmgpy/reactionTest.py \
  --no-cov -p no:cacheprovider -q -m "not functional and not database"
# full RMG-Py, two collectable halves
python -m pytest test/ --ignore=test/rmgpy/data/rmgTest.py --no-cov -p no:cacheprovider -q
python -m pytest test/rmgpy/data/rmgTest.py --no-cov -p no:cacheprovider -q
# i134 at tip, and at base (worktree of 71ae97bd5, rmgrc -> /home/alon/Code/RMG-database-plasma/input)
python -m pytest test/rmgpy/i134DuplicateElectronsTest.py --no-cov -p no:cacheprovider -q
# RMG-database (run from ../RMG-database-plasma)
PYTHONPATH=/home/alon/Code/RMG-Py-i221-saturation-fix python -m pytest test/ -p no:cacheprovider -q
# end to end
MPLCONFIGDIR=$TMPDIR python docs/i221-saturation/probes/round112_end_to_end_probe.py
```
