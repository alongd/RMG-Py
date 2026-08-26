# I-123 — does the merge backlog integrate?

**Verdict: NOT READY.** One blocking failure, named in §6.1: on the integrated pair, RMG's
Chemkin writer refuses every reaction carrying a `VoronovEIArrhenius` or `BadnellRRArrhenius`
rate law, so a lithium plasma model that the reactor accepts and evaluates correctly cannot be
saved. `rmg.py` dies inside `save_everything()` on the first model-save of the very deck the
backlog exists to enable.

Two integration branches were built, both from `origin/plasma`, ten merges total, **all
conflict-free**. Nothing was merged to `plasma`, nothing was pushed, no PR was opened.

| | |
|---|---|
| RMG-Py | `/home/alon/Code/RMG-Py-i123-integration`, branch `i123-integration`, tip `ae52de38a` |
| RMG-database | `/home/alon/Code/RMG-database-i123-integration`, branch `i123-integration-db`, tip `c7bd96292` |
| Build | `make build` — **exit 0**, on both the baseline and the union |
| RMG-Py unit suite | baseline **2596 passed / 0 failed**, union **2816 passed / 2 failed** |
| RMG-database suite | baseline **41 passed / 4 failed**, union **125 passed / 4 failed** (same 4) |
| Blocking failure | Chemkin writer refuses the plasma rate laws (§6.1) |

Evidence labels: **[R]** code (file:line) · **[D]** database (file, entry) · **[M]** measured
(command and output shown) · **[I]** inference with its basis · **UNKNOWN** where nothing was
established.

---

## 1. The two claims

### Claim 1 — the chain already collapsed: **CONFIRMED**, with one inaccuracy in its wording

```
$ for b in i093-isbalanced i097-library-electrons i101-fixtures i100-implicit-electrons \
           i108-electron-placement i113-placement-widening i116-ionisation-registry; do
    git merge-base --is-ancestor "$b" i119-rr-registry && echo "$b ANCESTOR"; done
i093-isbalanced              ANCESTOR  (tip aee9e510e)
i097-library-electrons       ANCESTOR  (tip 04bbb2881)
i101-fixtures                ANCESTOR  (tip 6cb64c70d)
i100-implicit-electrons      ANCESTOR  (tip 0a3b0ff3d)
i108-electron-placement      ANCESTOR  (tip 98b0e70b4)
i113-placement-widening      ANCESTOR  (tip dde602778)
i116-ionisation-registry     ANCESTOR  (tip 18abb2789)
```

**[M]** All seven are ancestors of `i119-rr-registry`. Those seven are one merge, not seven.

**[M]** The sub-claim "in a single linear history" is **wrong**, harmlessly. The range is 17
commits and two of them are merge commits:

```
$ git rev-list --merges origin/plasma..i119-rr-registry
42b8c29049bd8ac2e1424b4f58287a040f6b726d    # Merge the Marcus dGrxn fix
0c107ed379863b3e4b738f1f1a8adff56d50f432    # Merge walls A-D
```

It changes nothing about the merge order — noted so the next reader is not surprised by a
`--graph` that forks.

### Claim 2a — RMG-Py `i114-ionisation-declaration` is superseded by `i116-ionisation-registry`

**CONFIRMED on the thing that matters, REFUTED on the words "strict superset."** Skip `i114`.

The brief asked specifically whether the *functional* change is genuinely identical, not merely
similar, since a difference in the placement tuple or the scope rule would invert the answer. It
is identical. Both branches fork from `dde602778` and each adds exactly one commit.

**[M]** Filtering the `electron_placement.py` diff to non-comment lines leaves only the module
docstring's *scope* paragraph — prose, restating the same rule (a library label is a legal key
because `LibraryReaction` sets `family = library`). The declaration table itself is byte-identical:

```
$ git show i114-ionisation-declaration:rmgpy/electron_placement.py | grep -A4 '^FAMILY_ELECTRON_PLACEMENT'
FAMILY_ELECTRON_PLACEMENT = {                     # ← i114, line 152
    'Plasma_Electron_Attachment': (1, 0),
    'Cation_R_Recombination': (1, 0),
    'PlasmaElectronImpactIonization': (1, 2),
}
$ git show i116-ionisation-registry:rmgpy/electron_placement.py | grep -A4 '^FAMILY_ELECTRON_PLACEMENT'
FAMILY_ELECTRON_PLACEMENT = {                     # ← i116, line 149
    'Plasma_Electron_Attachment': (1, 0),
    'Cation_R_Recombination': (1, 0),
    'PlasmaElectronImpactIonization': (1, 2),
}
```

Same key, same `(1, 2)` tuple, same scope rule. **[R]** `rmgpy/electron_placement.py:152` (i114)
and `:149` (i116).

**The "strict superset" wording is false, and worth stating precisely.** `i114` carries one
assertion that `i116` deletes:

```python
# test/rmgpy/i113IonisationPlacementTest.py, on i114 only
assert 'PlasmaElectronImpactIonization' in FAMILY_ELECTRON_PLACEMENT
```

**[M]** No coverage is lost by that deletion: `i116` asserts the whole registry by equality four
lines earlier in the same class (`test/rmgpy/i113IonisationPlacementTest.py:402`), which is
strictly stronger — membership plus the value plus closure. And `i116` adds
`test/rmgpy/i116IonisationRegistryTest.py`, 650 lines that `i114` has no counterpart for.

**[M]** Merging `i114` after `i116` is not a no-op and not a prose revert: it conflicts.

```
$ git merge-tree --write-tree i116-ionisation-registry i114-ionisation-declaration; echo $?
1
CONFLICT (content): Merge conflict in docs/i108-electron-representation-remeasurement.md
CONFLICT (content): Merge conflict in rmgpy/electron_placement.py
CONFLICT (content): Merge conflict in test/rmgpy/electronPlacementTest.py
CONFLICT (content): Merge conflict in test/rmgpy/i108ElectronRepresentationMatrixTest.py
CONFLICT (content): Merge conflict in test/rmgpy/i113IonisationPlacementTest.py
```

Both branches rewrote the same docstrings differently. So skipping `i114` is right for a
stronger reason than the brief gave: not "it would only revert prose", but "it would cost five
hand-resolved conflicts to arrive at code that is already there."

`i116-ionisation-registry` arrives anyway as an ancestor of `i119-rr-registry` (Claim 1).

### Claim 2b — RMG-database `i102-quarantine-db` is superseded by `i111-sei-reclassification`

**CONFIRMED, on exactly the test the brief asked for.** Skip `i102-quarantine-db`.

The brief warned that a superset by line count is not a superset, and required confirmation that
no line present in the first version is absent from the second. Checked line by line:

```
$ git show i102-quarantine-db:input/kinetics/families/Cation_R_Recombination/quarantine.py > /tmp/q102.py
$ git show i111-sei-reclassification:input/kinetics/families/Cation_R_Recombination/quarantine.py > /tmp/q111.py
$ wc -l /tmp/q102.py /tmp/q111.py
  49 /tmp/q102.py
  94 /tmp/q111.py
$ diff -u /tmp/q102.py /tmp/q111.py
@@ -17,6 +17,40 @@   (34 added lines: classification, domain, isPlasmaFamily,
                       familySets, provenance, and the I-111 header)
@@ -46,4 +80,15 @@   (11 added lines: the I-111 note on where the gate exists)
```

**[M]** Every hunk is pure addition. Zero deleted lines, zero modified lines: `i111`'s version
contains `i102-quarantine-db`'s 49 lines verbatim and adds 45. **[D]**
`input/kinetics/families/Cation_R_Recombination/quarantine.py`. `i111` additionally reclassifies
`groups.py` and `recommended.py`, which `i102-quarantine-db` does not touch at all.

### A third collapse the brief did not claim

**[M]** RMG-database `i114-ionisation-declaration` is an **ancestor** of `i119-recombination`:

```
$ git merge-base --is-ancestor i114-ionisation-declaration i119-recombination && echo ANCESTOR
ANCESTOR
```

Database steps 4 and 5 are one chain, not two independent branches. Both merges were performed
anyway, in the stated order, so a failure could be attributed to the ionisation half separately
from the recombination half. The second merge brought only the recombination delta.

### The two "live work" exclusions currently exclude nothing

**[M]** Both excluded branches are at the identical commit as a branch being merged:

```
$ git rev-parse i119-rr-registry i121-transport-probe        # RMG-Py
92e2d423499504c643f563deb3e20b69797db186
92e2d423499504c643f563deb3e20b69797db186
$ git rev-parse i119-recombination i120-argon-recombination  # RMG-database
aabc3c622c85c1aac5e4b35966148d20e284e076
aabc3c622c85c1aac5e4b35966148d20e284e076
```

Neither has committed past its base yet. The exclusions cost nothing today and remain correct as
standing instructions; noted so nobody reads "excluded" as "its content is absent."

---

## 2. `origin/plasma` is not the plasma mainline

**[M]** The brief says both integration branches start from `origin/plasma`, and they do. But
`origin` is the personal fork `github.com/alongd/RMG-Py`, and its `plasma` ref is **ten commits
behind** the local `plasma` branch that `/home/alon/Code/RMG-Py-plasma` has checked out:

```
$ git rev-parse --short origin/plasma plasma
6d3c03823
a61dc1303
$ git log --oneline origin/plasma..plasma
a61dc1303 Merge the reverse-reconstruction refusal (I-098)
54efacb24 Refuse the reverse rate a reactor cannot reconstruct from Keq(Tgas)
fa3b58f05 Pin database.directory to the plasma database, relatively
42b8c2904 Merge the Marcus dGrxn fix
0c107ed37 Merge walls A-D: the corrected charged chemistry composes end to end
ce58b97f6 Resolve Marcus dGrxn from the reaction instead of the pressure
ab19481c8 Resolve electron placement for the reverse-generated cation reaction (I-088, Wall C)
0f5934dd2 Fix the electron sign on reverse-generated cation reactions (I-086, Wall D)
8857c3516 Add the M7-preflight reproduction and I-085 walls-A/B findings
c61504075 Clear the first two walls between a plasma input file and a running reactor
$ git log --oneline plasma..origin/plasma
                                            # empty: origin/plasma is strictly behind
```

**[I]** Consequence for this audit, and it is benign: those ten commits are ancestors of
`i115-preflight-deck` and `i112-marcus-work-terms`, so the union contains all of them anyway.
Verified: `git merge-base --is-ancestor plasma i115-preflight-deck` → yes. The baseline
number in §5 is therefore measured against a mainline ten commits staler than the real one — it
understates the baseline, and would have mattered had the union been red for a reason those ten
commits already carried. It is not (§6), so the delta stands as reported.

RMG-database has no such divergence: `origin/plasma` and `plasma` are both `fb3c13c60`.

---

## 3. The merge order actually used

Real merge commits, one per contributing branch, no squash and no rebase. Every branch was left
exactly as found.

### RMG-Py — `git log --oneline --graph --first-parent`

```
* ae52de38a Point the integration branch's rmgrc at the I-123 integration database
* f087cdf04 Merge i102-quarantine: hard-fail when quarantined database data reaches a reaction model
* 5e1686bcb Merge i112-marcus-work-terms: let the Marcus work terms reach the barrier they belong in
* 9f0eae40b Merge i115-preflight-deck: exclude the battery-SEI family from the plasma preflight decks
* c1f23e93e Merge i119-rr-registry: the charge-network chain, seven tickets in one branch
* 582fd6ab9 Merge i110-make-guard: refuse a bare `make` that installs into the shared conda environment
* 6d3c03823 Merge M4-B: the plasmaReactor input-file path with electron number-density conversion   ← origin/plasma
```

Full graph including second parents:

```
$ git log --oneline --graph origin/plasma..HEAD
*   ae52de38a Point the integration branch's rmgrc at the I-123 integration database
*   f087cdf04 Merge i102-quarantine
|\
| * 99602b099 docs: evidence for the Cation_R_Recombination Marcus quarantine
| * 541e6498f kinetics: hard-fail when quarantined database data reaches a reaction model
*   5e1686bcb Merge i112-marcus-work-terms
|\
| * ab3d84072 Let the Marcus work terms reach the barrier they belong in
*   9f0eae40b Merge i115-preflight-deck
|\
| * ceb9fd5f1 Exclude the battery-SEI family from the plasma preflight decks
| * a61dc1303 Merge the reverse-reconstruction refusal (I-098)     ← local plasma tip, see §2
| * ...
*   c1f23e93e Merge i119-rr-registry
|\
| * 92e2d4234 Let a kinetics library own radiative recombination
| * 18abb2789 Let a kinetics library own electron-impact ionisation   ← i116
| * dde602778 Let a family declare electrons on both sides            ← i113
| * ...                                                               ← i108, i100, i101, i097, i093
*   582fd6ab9 Merge i110-make-guard
|\
| * fa042c08d Document the shared-environment guard in CLAUDE.md
| * 5094dbea3 Refuse to install into the shared environment from a bare `make`
```

**Deviation from the brief, one, and it is the only non-merge commit on the branch:**
`ae52de38a` repoints `rmgrc`. The `rmgrc` that arrives with `i119-rr-registry` (commit
`fa3b58f05`) is a *tracked* file at the repo root pinning `database.directory` to
`../RMG-database-plasma/input`. Left as merged, the RMG-Py union would have been tested against
the **plasma** database worktree, which carries neither `PlasmaElectronImpactIonization` nor
`PlasmaRadiativeRecombination` — the two libraries this whole audit exists to exercise — and the
suite would have reported a green that meant nothing. The commit points it at
`../RMG-database-i123-integration/input` instead and says in its own message that it must be
reverted before any merge to `plasma`.

**[M]** Commit counts differed slightly from the brief in one place: `i102-quarantine` is 10
commits above `origin/plasma`, not 11. Every other count matched.

### RMG-database — `git log --oneline --graph`

```
*   c7bd96292 Merge i119-recombination: give radiative recombination a library owner
|\
| * aabc3c622 docs: report on the electron loss channels
| * 9ee5042d1 test: pin the recombination path from load to reactor acceptance
| * 9325c2d27 kinetics: give radiative recombination a library owner
* | 1fb224371 Merge i114-ionisation-declaration: give electron-impact ionisation a library owner
|\|
| * fcd87ce6e docs: record that the placement declaration landed
| * 277a0e333 test: pin the ionisation path from load to reactor acceptance
| * 6862c21c5 kinetics: give electron-impact ionisation a library owner
| * bdbbf2864 test: pin database.directory for the suites in this repository
*   ac944618e Merge i104-alkali-plasma: inventory the alkali-plasma cation source and loss network
|\
| * fc7bc101b docs: inventory the alkali-plasma cation source and loss network
*   79320d363 Merge i103-electrochem-provenance: diagnose the Cation_R_Recombination Marcus provenance
|\
| * 069ea3d0b docs: diagnose the Cation_R_Recombination Marcus provenance
* 3b61056c1 Merge i111-sei-reclassification: reclassify Cation_R_Recombination as a legacy SEI family
* 1fcbb80fa docs: report on the SEI reclassification
* aaa97074c test: pin the SEI reclassification and its plasma exclusion
* 03b96a974 kinetics: reclassify Cation_R_Recombination as a legacy SEI family
* 7d4f47d3a kinetics: quarantine the Cation_R_Recombination Marcus data
```

The `|\|` at `1fb224371` is the third collapse: `i114-ionisation-declaration` is an ancestor of
`i119-recombination`, so the two merges share a side of the history.

**No conflicts, in either repository, at any of the ten steps.** The brief said a conflict would
itself be a finding; there were none to report.

---

## 4. Build — stated before any test result

**[M]** `make build`, on the integrated RMG-Py tree, after all five merges:

```
$ make build
...
copying build/lib.linux-x86_64-cpython-39/rmgpy/solver/plasma.cpython-39-x86_64-linux-gnu.so -> rmgpy/solver
MAKE_BUILD_UNION_EXIT=0
```

**Exit 0.** 104 `.so` files in the worktree. The seventeen Cython-touching commits compiled
clean; no warnings escalated, no module failed to build. The same command on the baseline also
exited 0, so the build result is not a delta.

**[M]** The guard `i110-make-guard` adds works, checked immediately after step 1 as the brief
asked:

```
$ make
Refusing to modify the shared RMG environment.
Use `make build` for an in-place worktree build.
Editable installation requires an explicit maintenance procedure:
    make unsafe-install-shared-env CONFIRM_SHARED_ENV_MUTATION=yes
make: *** [Makefile:72: guard] Error 1
BARE_MAKE_EXIT=2
```

Bare `make` fails by design. Nothing repointed the shared conda environment at any point in this
audit.

**[M]** The environment resolved to the integration worktrees throughout:

```
$ python -c "import rmgpy, rmgpy.molecule.molecule as m; print(rmgpy.__file__); print(m.__file__)"
/home/alon/Code/RMG-Py-i123-integration/rmgpy/__init__.py
/home/alon/Code/RMG-Py-i123-integration/rmgpy/molecule/molecule.cpython-39-x86_64-linux-gnu.so
$ python -c "from rmgpy import settings; print(settings['database.directory'])"
/home/alon/Code/RMG-database-i123-integration/input
```

Never `/home/alon/Code/RMG-Py`. One database-repo run was launched without `PYTHONPATH`, resolved
to the shared checkout, and failed with `ModuleNotFoundError: No module named
'rmgpy.electron_balance'` — caught immediately by that error and re-run correctly; only the
correct run is reported below.

---

## 5. Test counts, both sides of the delta

Same command each time. RMG-Py suite run from the RMG-Py worktree root; database suite run from
the database worktree root with `PYTHONPATH` set to the RMG-Py integration worktree.

### RMG-Py — `python -m pytest -m "not functional and not database"`

| | passed | failed | skipped | deselected |
|---|---|---|---|---|
| `origin/plasma` baseline | **2596** | **0** | 49 | 68 |
| integrated pair, merges only | **2816** | **2** | 50 | 130 |
| integrated pair, final (this branch's tip) | **2818** | **2** | 49 | 131 |

```
baseline:      === 2596 passed, 49 skipped, 68 deselected, 54 warnings in 261.80s (0:04:21) ===
merges only:   = 2 failed, 2816 passed, 50 skipped, 130 deselected, 54 warnings in 243.61s (0:04:03) =
final:         = 2 failed, 2818 passed, 49 skipped, 131 deselected, 54 warnings in 347.85s (0:05:47) =
```

**+220 passing tests** — the suites the merged branches bring. **+2 failures**, both named and
attributed in §6, identical in both union runs:

```
FAILED test/rmgpy/electronPlacementTest.py::TestElectronPlacementResolver::test_declaration_registry_is_explicit_and_closed
FAILED test/rmgpy/preflightDeckFamilyExclusionTest.py::PlasmaDeckFamilyExclusionTest::test_plasma_deck_does_not_declare_family[docs/i102-quarantine/input.py]
```

The final row is the branch as it stands, and the +2/+1 over the merges-only row is this audit's
own deck (§7) being discovered by `i115-preflight-deck`'s repository-wide sweep: two more
parametrisations, both passing, one of them database-marked and therefore deselected. It changes
no failure. (`pytest` exits 0 on the union because the default `addopts` do not fail the run on a
coverage threshold; the `2 failed` line is the authority, not the exit code.)

### RMG-database — `python -m pytest test/ -q`

| | passed | failed | exit |
|---|---|---|---|
| baseline pair | **41** | **4** | 1 |
| integrated pair | **125** | **4** | 1 |

**+84 passing tests**, all from the ionisation, recombination and SEI suites. **The same four
failures on both sides**, so none is union-only:

```
FAILED test/test_plasma_electron_attachment.py::test_trained_species_resolve_to_their_own_rate_through_their_own_node[O2]
FAILED test/test_plasma_electron_attachment.py::test_trained_species_resolve_to_their_own_rate_through_their_own_node[OH]
FAILED test/test_plasma_electron_attachment.py::test_trained_species_resolve_to_their_own_rate_through_their_own_node[O]
FAILED test/test_plasma_electron_attachment.py::test_o2_rate_comes_from_the_training_set_not_a_library
```

All four assert `source == 'rate rules'` and get
`<KineticsDepository "Plasma_Electron_Attachment/training">`. **[M]** Pre-existing on the
baseline pair; **out of scope for this ticket** and not diagnosed here beyond establishing that
the merges did not cause them.

**[M]** A note on how the database baseline had to be taken. Run as-is at `origin/plasma`, the
database suite produces **45 errors and 0 passes**, because `test/conftest.py` — which pins
`database.directory` — arrives with `i114-ionisation-declaration` and does not exist on
`origin/plasma`. The baseline above was taken with an equivalent `rmgrc` at the database
worktree root so the two runs measure the same thing. Without that, the "delta" would have been
45 errors → 0, which is a measurement artefact and not a result.

---

## 6. Failures that exist only on the union

Three, in descending order of consequence. One is blocking; one is a genuine two-branch
interaction; one turns out not to be a union failure at all.

### 6.1 BLOCKING — the Chemkin writer refuses the plasma rate laws

**Not a pytest failure.** The unit suites are blind to it. It appears the moment a real model
containing one of these reactions is saved.

**[M]** Running the lithium deck (§7) through `rmg.py` on the integrated tree:

```
$ python rmg.py docs/i123-integration/input.py
...
Saving current model core and edge to Chemkin file...
Traceback (most recent call last):
  File ".../rmgpy/rmg/main.py", line 961, in execute
    self.save_everything()
  File ".../rmgpy/rmg/main.py", line 2254, in save_everything
    self.notify()
  File "rmgpy/chemkin.pyx", line 1829, in rmgpy.chemkin.write_reaction_string
    check_electron_reactant_order(reaction, reactants, reaction_string)
  File ".../rmgpy/electron_balance.py", line 226, in check_electron_reactant_order
    raise MechanismWriterError(
rmgpy.exceptions.MechanismWriterError: Reaction Li(2) => [Lip](3) has kinetics
VoronovEIArrhenius whose rate is proportional to the electron density, but the exported
equation "Li(2)=>[Lip](3)+e-(1)" has no electron among its reactants, so a solver would
evaluate it at the wrong reaction order. Put the consumed electron in reaction.reactants
and count only surplus produced electrons in reaction.electrons.
```

**The two designs that collide, both correct in isolation:**

1. **[R]** `rmgpy/electron_balance.py:203-232`, `check_electron_reactant_order` — a write-time
   guard: a reaction whose rate law is proportional to electron density must carry the consumed
   electron among `reaction.reactants`, or the exported file is silently wrong by a factor of the
   electron density. Called from `rmgpy/chemkin.pyx:1829`.

2. **[D]** `input/kinetics/libraries/PlasmaElectronImpactIonization/reactions.py` — the canonical
   database form is `[Li] => [Lip]` with `electrons = +1` and **no explicit electron
   participant**, deliberately: *"the entry below is the canonical database form and must not be
   written with an explicit electron participant (the resolver refuses double representation)."*
   The electron is restored only into a non-mutating reactor-facing *view* produced by
   `resolve_electron_placement`. Same for `PlasmaRadiativeRecombination`.

**[M]** The join is missing. The only production call site of the resolver is inside the reactor:

```
$ grep -rn "resolve_electron_placement" --include=*.py --include=*.pyx rmgpy/ | grep -v test
rmgpy/solver/plasma.pyx:236   core_reactions = self._resolve_electron_placements(...)
rmgpy/solver/plasma.pyx:237   edge_reactions = self._resolve_electron_placements(...)
rmgpy/solver/plasma.pyx:334       electron_placement.resolve_electron_placement(rxn, species_list))
rmgpy/electron_placement.py:267   def resolve_electron_placement(reaction, species_list):
```

`PlasmaReactor.initialize_model` resolves placement for itself. The Chemkin writer never does,
and receives the canonical reaction. So the reactor accepts what the writer refuses, and RMG
calls the writer on every model save.

**Attribution — and it is worse than a two-branch interaction.** **[M]** The guard is *not*
introduced by any merged branch:

```
$ git log --oneline --all -S 'check_electron_reactant_order' -- rmgpy/electron_balance.py
a9c81e6de0 Export plasma reactions instead of dropping them, and carry the electron
$ git diff --stat origin/plasma HEAD -- rmgpy/electron_balance.py rmgpy/chemkin.pyx
                                            # empty: untouched by all five merges
```

It is already on `origin/plasma`. What the union adds is data that trips it. Minimal reproducer,
with **both** of `i119-rr-registry`'s registry entries deleted at runtime, to prove the RMG-Py
half is not required:

```
$ python  # registry stripped to {'Plasma_Electron_Attachment', 'Cation_R_Recombination'}
canonical reaction : [Li] => [Lip] | reactants: ['[Li]']
RESULT: MechanismWriterError: ... has electrons=1 but the mechanism defines no electron
species, so the electron cannot be written into the exported equation.
```

**[I]** The refusal survives with the RMG-Py registry entries removed — it is triggered by the
library *data* alone against a guard that is already on the mainline. The interacting pair is
therefore **RMG-database `i114-ionisation-declaration` (equivalently its descendant
`i119-recombination`) × RMG-Py `plasma` as it already stands**, not two of the branches under
audit. **Merging the RMG-database ionisation branch alone into `plasma` would break Chemkin
export for any plasma model that uses it.** RMG-Py `i119-rr-registry` is what makes the reaction
reactor-*legal*, so it is what makes the deck get far enough to hit the writer; it is not what
makes the writer refuse.

**Not fixed here.** The brief forbids changing `rmgpy/chemkin.pyx` or `rmgpy/electron_balance.py`
beyond what the merges bring, and forbids fixing a failing test by changing it. The red is
reported. **[I]** The design question a fix has to answer — and it is a real one, not a typo —
is *which* representation the exported mechanism should carry, since the guard and the library
were each written against a different answer. That is a ticket, not an edit.

### 6.2 GENUINE union-only failure — a deck sweep meets a new deck

**Two failures, from one interaction**, and only one of them is visible to the unit suite:

```
FAILED ...::test_plasma_deck_does_not_declare_family[docs/i102-quarantine/input.py]
E  AssertionError: plasma deck docs/i102-quarantine/input.py declares Cation_R_Recombination,
   a lithium-ion-battery SEI family excluded from every plasma configuration.
   Declared families: ['Cation_R_Recombination', 'Cation_Li_Abstraction', ...]

FAILED ...::test_plasma_deck_generates_no_reaction_from_family[docs/i102-quarantine/input.py]
```

**[M]** The second carries `@pytest.mark.database`
(`test/rmgpy/preflightDeckFamilyExclusionTest.py:230`) and is therefore **deselected** from
`-m "not functional and not database"`. It appears only when the file is run directly. The union
unit-suite count in §5 undercounts this interaction by one, and §10.4 is where that generalises.

**Interacting pair: `i115-preflight-deck` × `i102-quarantine`.** Each is green alone.

- **[R]** `i115-preflight-deck` adds `test/rmgpy/preflightDeckFamilyExclusionTest.py`, which
  **discovers plasma decks by AST-sweeping the whole repository** for files calling
  `plasmaReactor` (`_discover_decks()`, line 165) and asserts that no discovered deck declares
  `Cation_R_Recombination`. On its own branch it finds three decks, all under `docs/m7-preflight/`,
  and all three pass.
- **[R]** `i102-quarantine` adds `docs/i102-quarantine/input.py`, a plasma deck that declares
  `Cation_R_Recombination` — legitimately, since its whole purpose is to demonstrate the
  quarantine gate firing on that family.

Neither branch can see the other. The sweep only finds the fourth deck once both are merged. This
is exactly the class of failure the brief predicted would be worth the most, and the only one of
the three that is created by merging two of the audited branches together.

**[M]** The other four parametrisations still pass (`docs/m7-preflight/input.py`,
`input_placement.py`, `input_secondary.py`, and the SEI negative control), so the sweep is working
and the failure is specific.

**[I]** This one is cheap to resolve and the choice is a judgement call, not a defect: either
`docs/i102-quarantine/input.py` is exempted by name as a deliberate gate-demonstration deck, or
the sweep is narrowed to decks intended for execution. Left unfixed, as instructed.

### 6.3 NOT a union failure — `i119-rr-registry` is red on its own branch

```
FAILED test/rmgpy/electronPlacementTest.py::TestElectronPlacementResolver
       ::test_declaration_registry_is_explicit_and_closed
E  assert {'Plasma_Electron_Attachment': (1,0), 'Cation_R_Recombination': (1,0),
          'PlasmaRadiativeRecombination': (1,0), 'PlasmaElectronImpactIonization': (1,2)}
       == {'Plasma_Electron_Attachment': (1,0), 'Cation_R_Recombination': (1,0),
          'PlasmaElectronImpactIonization': (1,2)}
E  Left contains 1 more item: {'PlasmaRadiativeRecombination': (1, 0)}
```

The brief's opening premise is that all eleven branches are "committed, green on their own
branches." **[M]** For `i119-rr-registry` that is false, and the merge imported the failure
rather than causing it.

Its tip commit `92e2d4234` ("Let a kinetics library own radiative recombination") adds a fourth
registry entry. There are **three** closed-registry assertions in the suite. It updated two:

```
$ git show 92e2d4234 --stat
 rmgpy/electron_placement.py                        | 32 +++++++++++++++++++---
 test/rmgpy/i108ElectronRepresentationMatrixTest.py | 17 +++++++++---
 test/rmgpy/i113IonisationPlacementTest.py          | 20 ++++++++++----
```

`test/rmgpy/electronPlacementTest.py` is not in that list, and still asserts three entries.

**[M]** Measured directly against `i119-rr-registry`'s own two files, parsed out of the branch
with no merge involved:

```
i119-rr-registry rmgpy/electron_placement.py  FAMILY_ELECTRON_PLACEMENT =
    {'Plasma_Electron_Attachment': (1,0), 'Cation_R_Recombination': (1,0),
     'PlasmaRadiativeRecombination': (1,0), 'PlasmaElectronImpactIonization': (1,2)}
i119-rr-registry test/rmgpy/electronPlacementTest.py asserts equality against =
    {'Plasma_Electron_Attachment': (1,0), 'Cation_R_Recombination': (1,0),
     'PlasmaElectronImpactIonization': (1,2)}
EQUAL? False
registry has, that the assertion lacks: {'PlasmaRadiativeRecombination'}
```

**[M]** And both files entered the union byte-identically — the merge changed neither:

```
rmgpy/electron_placement.py         i119=249298e38c union=249298e38c IDENTICAL
test/rmgpy/electronPlacementTest.py i119=617d58c389 union=617d58c389 IDENTICAL
```

`rmgpy/electron_placement.py` is not cythonized, so no stale `.so` is involved. The failure is a
one-line omission on `i119-rr-registry`: add `'PlasmaRadiativeRecombination': (1, 0)` to the
assertion at `test/rmgpy/electronPlacementTest.py:351`. **Not fixed here** — it belongs on
`i119-rr-registry`, and this branch merges rather than edits.

**[I]** The wider point: "green on its own branch" was assumed, not measured, for eleven
branches. It has now been measured for one, and was wrong. §10 names this as the largest thing
this audit did not reach.

---

## 7. The lithium charge network, end to end on the integrated tree

Two independent demonstrations. The probe drives the loader, the resolver and `PlasmaReactor`
directly, one stage at a time; the deck drives the same chemistry through the user-facing
`rmg.py` path. Both are committed: `docs/i123-integration/probe_lithium_charge_network.py` and
`docs/i123-integration/input.py`.

**Why not `docs/m7-preflight/input.py`.** **[R]** That deck's own header, after
`i115-preflight-deck` removed `Cation_R_Recombination` from it, records that it now generates no
reactions at all: the reverse of that family's template was its only route from the neutral feed
to Li⁺, and *"restoring it needs a genuine plasma ionisation route (e.g. electron-impact
ionisation), not this family."* That route is precisely what the union supplies. The new deck is
that deck's stated successor.

### The probe — every stage separately

**[M]** `python docs/i123-integration/probe_lithium_charge_network.py` → exit 0, all 22 checks
pass. Full log: `docs/i123-integration/evidence/probe-lithium-charge-network.log`.

**Stage 1 — LOAD.** Both libraries load out of the integration database.

```
ionisation    : [Li] => [Lip]   electrons=1   kinetics=VoronovEIArrhenius
recombination : [Lip] => [Li]   electrons=-1  kinetics=BadnellRRArrhenius
[PASS] ionisation declares a net electron GAIN -- 1
[PASS] recombination declares a net electron LOSS -- -1
[PASS] the neutral feed and the cation are the same two species in both channels
       -- source produces exactly what the sink consumes
```

The last check is the integration claim in one line: the species the source produces is
isomorphic to the one the sink consumes, so the two libraries describe one network and not two.

**Stage 2 — BALANCE.** `Reaction.is_balanced` — which only became a real charge check on
`i093-isbalanced` — closes both gaps:

```
[PASS] ionisation is_balanced          (0 -> +1, closed by electrons=+1)
[PASS] recombination is_balanced       (+1 -> 0, closed by electrons=-1)
```

**Stage 3 — PLACEMENT.** The resolver turns the canonical form into the reactor-facing one, from
the declarations that `i116` (ionisation) and `i119` (recombination) added:

```
FAMILY_ELECTRON_PLACEMENT = {'Plasma_Electron_Attachment': (1, 0),
                             'Cation_R_Recombination': (1, 0),
                             'PlasmaRadiativeRecombination': (1, 0),
                             'PlasmaElectronImpactIonization': (1, 2)}
ionisation view    : [Li] + e => [Lip] + e + e  -> reactant/product electrons (1, 2), order 2
recombination view : [Lip] + e => [Li]          -> reactant/product electrons (1, 0), order 2
[PASS] both resolved views are second order in the reactor
[PASS] resolution did not mutate the canonical reactions
```

Both views are second order, which is what makes them dimensionally consistent with the
`cm^3/(molecule*s)` A-factors, and the canonical reactions are unchanged — the non-mutation
contract holds.

**Stage 4 — REACTOR.** One `PlasmaReactor` holds the source *and* the sink, resolves placement
itself, and evaluates each at its own rate law:

```
[PASS] the reactor accepted both reactions -- 2
[PASS] the reactor located the electron -- 2
kf[source] = 5.031584e+07   expected (Voronov) = 5.031584e+07
kf[sink]   = 1.468188e+05   expected (Badnell) = 1.468188e+05
[PASS] the cation has a non-zero loss channel
```

Not merely accepted — evaluated at the right numbers, to relative 1e-9, at Te = 10 000 K. A view
resolved at the wrong reaction order would be silent, and the `kf` comparison is what makes it
visible.

**This is the first time all of it has run with every piece present at once, and it works.**

### The deck — the same chemistry through `rmg.py`

**[M]** `python rmg.py docs/i123-integration/input.py`, neutral atomic Li feed at 0.05 in He, no
cation seeded, Te = 11 604.5 K (~1 eV), both libraries in `reactionLibraries`, and two neutral
lithium families declared — neither of which produces a cation, so the charge network is
library-only.

**A trap found while writing this deck, recorded because the next deck will hit it.** The obvious
way to say "no families" is `kineticsFamilies=[]`, and it does the opposite: an empty list is
falsy, and `KineticsDatabase.load_families` reads a falsy `families` as *load the default set*, so
the deck silently gets **every** family — including `Cation_R_Recombination`, the quarantined
battery-SEI family that `i111` and `i102` exist to keep out of plasma decks. **[M]** Measured, not
inferred:

```
$ pytest "...::test_plasma_deck_generates_no_reaction_from_family[docs/i123-integration/input.py]"
E  AssertionError: plasma deck docs/i123-integration/input.py resolves to a family set that
   loads Cation_R_Recombination. Declared: []
1 failed in 189.22s
```

**[I]** `kineticsFamilies='none'` loads nothing, but the deck-exclusion sweep asserts
`isinstance(families, list)`, so `'none'` fails a different assertion. There is currently no
spelling of "a plasma deck with no families" that satisfies both. Out of scope here; the deck
names an explicit two-family list instead, and **[M]** now passes both sweep tests:

```
$ pytest test/rmgpy/preflightDeckFamilyExclusionTest.py -k i123
test_plasma_deck_does_not_declare_family[docs/i123-integration/input.py] PASSED
test_plasma_deck_generates_no_reaction_from_family[docs/i123-integration/input.py] PASSED
================= 2 passed, 13 deselected in 172.64s (0:02:52) =================
```

What the deck proves before it dies:

```
Global RMG Settings:
   database.directory   = /home/alon/Code/RMG-database-i123-integration/input (from rmgrc)
The current git HEAD for RMG-Py is:       ae52de38ae1fd552efe3af9108bd52f20f4d28f1
The current git HEAD for RMG-database is: c7bd96292ef92d58a0e2c7a68a6e5ebc93188a26
...
    The model core has 6 species and 0 reactions
    The model edge has 1 species and 1 reactions
```

Six core species — `Li`, `He`, `e-`, `[Lip]` and their companions — and the edge holds exactly the
one ionisation reaction, with no family-generated noise.

The input file parses, both libraries load, the electron pseudo-species is declared and accepted,
`[Lip]` is constructed from the neutral feed and gets thermochemistry, and the ionisation reaction
enters the model — the reactor path all the way to the first enlargement.

Then it dies in `save_everything()` on §6.1. **[M]** The failure is on the *export* path, after
the reactor path has succeeded; the probe (which never exports) completes all four stages. Both
statements are true and they are about different halves of the system.

Evidence: `docs/i123-integration/evidence/deck-run-chemkin-writer-refusal.{stdout.tail,stderr}.log`.

---

## 8. Negative control — an ordinary non-plasma mechanism is unaffected

**[M]** `python rmg.py examples/rmg/minimal/input.py` on the integrated tree:

```
MODEL GENERATION COMPLETED
The final model core has 26 species and 66 reactions
The final model edge has 146 species and 332 reactions
RMG execution terminated at Wed Aug 26 04:50:41 2026
NEGCTL_EXIT=0
```

Exit 0, 27 seconds, zero tracebacks in `stderr.log`, `chem.inp` and `chem_annotated.inp` both
written. Nothing the union adds — not the quarantine gate `i102-quarantine` installs, not the
Marcus work terms, not the electron-placement machinery — perturbs an ordinary gas-phase run.

The quarantine gate is the one worth calling out: `i102-quarantine` makes quarantined database
data a hard failure when it reaches a reaction model, and an over-broad gate would have shown up
here as a refusal. It did not fire. Evidence:
`docs/i123-integration/evidence/negative-control.log`.

---

## 9. Verdict

**NOT READY**, on §6.1 alone.

| Finding | Blocking? | Owner |
|---|---|---|
| §6.1 Chemkin writer refuses `VoronovEIArrhenius` / `BadnellRRArrhenius` | **yes** | needs a ticket: which representation does an exported plasma mechanism carry? |
| §6.2 deck sweep × quarantine deck | no | one-line exemption or a narrowed sweep; a judgement call |
| §6.3 `i119-rr-registry` red on its own branch | no | one line on `i119-rr-registry`, not here |
| §5 four pre-existing `Plasma_Electron_Attachment` failures | no | pre-existing, out of scope |

Everything else the audit could check is green: ten conflict-free merges, a clean build, +220
and +84 passing tests, the lithium charge network working end to end through the reactor, and an
untouched negative control.

**On whether the gate has been protecting the mainline or deferring the cost — it deferred it.**
The merges themselves were free: no conflicts, no rebases, nothing to hand-resolve. What was
deferred is one real defect that no branch could have found alone, because it needs RMG-database
data and RMG-Py code in the same process, and because it lives on the export path that neither
repository's unit suite reaches. That defect is now visible, three lines of traceback deep,
instead of being discovered by whoever merged first.

**§6.1 also changes what "merge one branch at a time" is worth here.** It is not gated on the
full backlog: merging RMG-database `i114-ionisation-declaration` on its own breaks Chemkin export
for any plasma model that uses the library. That is the one ordering fact in this report that
should reach whoever merges next.

---

## 10. What this audit did not reach, and what would settle it

1. **Whether the other ten branches are green on their own branches.** Only `i119-rr-registry` was
   checked, because its failure surfaced in the union and demanded attribution — and it was
   **not** green (§6.3). The premise is untested for the rest. Settled by: checking out each
   branch and running `make build && make test` against a correctly pinned database, ten times.
   That is a second ticket's worth of compute, and §6.3 is reason enough to spend it.

2. **Whether §6.1 also breaks the Cantera and RMS export paths.** Only the Chemkin writer was
   reached, because it is the first writer `save_everything()` calls and the run dies there.
   `rmgpy/yaml_cantera2.py` and `rmgpy/yaml_rms.py` are both touched by `i112-marcus-work-terms`
   and both export reactions. **UNKNOWN.** Settled by: stubbing the Chemkin writer out of the
   observer list and re-running the deck, or unit-testing each writer against a canonical
   `VoronovEIArrhenius` reaction.

3. **Whether the deck's chemistry is right past the first enlargement.** The run never reached a
   converged model, so no rate, no time profile and no electron-density trajectory was checked.
   The probe pins `kf` at one temperature for each channel; nothing pins the *simulation*.
   **UNKNOWN.** Settled by: §6.1 being fixed, then running the deck to termination.

4. **Functional and database-marked tests — and this is not a formality.** The brief asked for
   the unit suite and that is what was run (`-m "not functional and not database"`); 130 tests
   were deselected on the union against 68 on the baseline, so the merges nearly doubled the
   deselected set and 62 newly-added tests were never executed. **[M]** One of them is already
   known to be red: `test_plasma_deck_generates_no_reaction_from_family[docs/i102-quarantine/input.py]`
   (§6.2), found only because that file was run directly. **[I]** One deselected failure found
   by accident in a set of 62 is a poor argument that the rest are green. **UNKNOWN.** Settled
   by: `make test-all` on the integrated pair. Note the cost — that one file alone takes 13
   minutes, because each parametrisation loads the full kinetics database.

5. **The four pre-existing `Plasma_Electron_Attachment` failures.** Established as pre-existing
   and not diagnosed. They assert `source == 'rate rules'` and receive a `KineticsDepository`, so
   the shape suggests a rate-rule-versus-training-depository precedence question rather than a
   chemistry error, but that is a guess and is labelled as one. **UNKNOWN.**

6. **Whether `origin/plasma` or local `plasma` is the intended merge target.** §2 establishes the
   ten-commit divergence as fact; which one is *supposed* to be the mainline is a question for
   whoever owns the branch, not a thing this audit can measure. The stale ref is on a personal
   fork, which suggests the local branch is authoritative and the fork simply was not pushed —
   but that is inference. **UNKNOWN.** Settled by: the branch owner saying so, or a push.

---

## Reproducing this

```bash
conda activate rmg_env
cd /home/alon/Code/RMG-Py-i123-integration
export PYTHONPATH=/home/alon/Code/RMG-Py-i123-integration:$PYTHONPATH
python -c "import rmgpy; print(rmgpy.__file__)"          # must be under this worktree
python -c "from rmgpy import settings; print(settings['database.directory'])"
                                                          # must be RMG-database-i123-integration/input

make build                                                # never bare `make`
python -m pytest -m "not functional and not database"     # §5, RMG-Py
python docs/i123-integration/probe_lithium_charge_network.py   # §7, all four stages
python rmg.py docs/i123-integration/input.py              # §6.1, dies in save_everything()
python rmg.py examples/rmg/minimal/input.py               # §8, negative control

cd /home/alon/Code/RMG-database-i123-integration
python -m pytest test/ -q                                 # §5, RMG-database
```

Long-running commands in this audit were run as
`<command> > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)`; the §6.1 traceback is only in
`stderr.log`, not in `RMG.log`, which is why it is captured.
