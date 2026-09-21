# I-244 — the Chemkin writer marked two reactions of different stoichiometry as duplicates

Worktree `/home/alon/Code/RMG-Py-i244-chemkin-duplicate-electron-aware`, branch
`i244-chemkin-duplicate-electron-aware`, cut from `plasma` at `311818121`.
Database `/home/alon/Code/RMG-database-plasma`, read-only and unmodified.
All runs captured to `logs/` with both streams.

---

## 1. What changed

One function, `mark_duplicate_reaction` in `rmgpy/chemkin.pyx`. The branch that **marks** a pair
now consults the per-side electron placement; the two branches that **un**-mark a pair are
untouched. Nothing else in the tree changed.

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

For every reaction outside the plasma families and libraries both counts are `(0, 0)`, so the
marking verdict is unchanged bit for bit. §5 measures that rather than arguing it.

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
