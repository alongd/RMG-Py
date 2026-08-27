# I-123, third pass — the union against the right baseline

**Verdict: NOT READY.** One blocking failure, §8: a plasma mechanism that RMG reads back and
writes out again carries every rate constant **exactly 10⁶ too small**, silently, under a units
header the numbers are not in. The reader repair merged in this pass is what makes a pre-existing
writer defect reachable, and the shipped `scripts/mergeModels.py` reaches it end to end with exit
status 0.

Everything the pass was asked to prove otherwise holds. Both previously-blocking failures are
genuinely fixed; the lithium charge network runs end to end and both channels reach the generated
model; the negative control is byte-identical; the merges are conflict-free.

| | |
|---|---|
| RMG-Py | `/home/alon/Code/RMG-Py-i123c-audit`, branch `i123c-audit` |
| RMG-database | `/home/alon/Code/RMG-database-i123c-audit`, branch `i123c-audit-db` |
| Authoritative base | **local** `plasma` = `a61dc1303` (RMG-Py) · local `plasma` = `fb3c13c60` (RMG-database) |
| Stale remote | `origin/plasma` = `6d3c03823` (RMG-Py), **10 commits behind local** · `fb3c13c60` (RMG-database, **identical**) |
| Build | `make build` exit **0**; bare `make` **refused** on the union |
| RMG-Py unit suite | base (local `plasma`) **2623 passed / 0 failed** · union **3020 passed / 2 failed** (both pre-existing, §6) |
| Blocking failure | §8 — silent 10⁶ rate corruption on Chemkin re-export of a reloaded plasma mechanism |
| Merged to a shared branch | **nothing**. Pushed: **nothing**. PR: **none**. |

Evidence labels: **[R]** code (file:line) · **[D]** database (file, entry) · **[M]** measured
(command and output shown) · **[I]** inference with its basis · **UNKNOWN** where nothing was
established. Raw output is under `docs/i123c-audit/evidence/`.

Everything below was run with `PYTHONPATH=/home/alon/Code/RMG-Py-i123c-audit`, against
`database.directory = /home/alon/Code/RMG-database-i123c-audit/input`, **printed at the head of
every probe** rather than trusted from configuration [M]:

```
rmgpy              = /home/alon/Code/RMG-Py-i123c-audit
compiled module    = /home/alon/Code/RMG-Py-i123c-audit/rmgpy/molecule/molecule.cpython-39-x86_64-linux-gnu.so
database.directory = /home/alon/Code/RMG-database-i123c-audit/input
resolved realpath  = /home/alon/Code/RMG-database-i123c-audit/input
exists             = True
```

---

## 1. The baseline correction — and what it actually was

The brief's premise: both previous audits measured against the **remote** shared branch, the
remote is stale by ten commits, so their baselines were wrong and any of the ten could collide
with the union. **Half of that is right, and the half that is wrong is the half that mattered.**

### 1.1 The ten commits, shown

```
$ git rev-parse --short plasma origin/plasma
a61dc1303
6d3c03823

$ git rev-list --count origin/plasma..plasma      # local ahead of remote
10
$ git rev-list --count plasma..origin/plasma      # remote ahead of local
0
```

**[M]** Ten commits, oldest first:

```
$ git log --reverse --format='%h %ad %s' --date=short origin/plasma..plasma
c61504075 2026-08-23 Clear the first two walls between a plasma input file and a running reactor
8857c3516 2026-08-23 Add the M7-preflight reproduction and I-085 walls-A/B findings
0f5934dd2 2026-08-23 Fix the electron sign on reverse-generated cation reactions (I-086, Wall D)
ab19481c8 2026-08-23 Resolve electron placement for the reverse-generated cation reaction (I-088, Wall C)
ce58b97f6 2026-08-24 Resolve Marcus dGrxn from the reaction instead of the pressure
0c107ed37 2026-08-24 Merge walls A-D: the corrected charged chemistry composes end to end
42b8c2904 2026-08-24 Merge the Marcus dGrxn fix
fa3b58f05 2026-08-24 Pin database.directory to the plasma database, relatively
54efacb24 2026-08-25 Refuse the reverse rate a reactor cannot reconstruct from Keq(Tgas)
a61dc1303 2026-08-25 Merge the reverse-reconstruction refusal (I-098)
```

**What they change.** Four substantive lines of work plus their merges and one configuration
commit: the M7-preflight walls A–D (electron sign on reverse-generated cations, electron placement
for the same, electron thermo), the Marcus `dGrxn` resolution, the reverse-reconstruction refusal
guard, and `fa3b58f05` which **created the tracked `rmgrc`** this lineage now depends on. 43 files,
+4199/−118.

### 1.2 Do any of them touch what the union touches? **Yes — ten files.**

```
$ comm -12 <(git diff --name-only origin/plasma...plasma | sort) \
           <(git diff --name-only plasma..HEAD | sort)
docs/m7-preflight/input_placement.py
docs/m7-preflight/input.py
docs/m7-preflight/input_secondary.py
rmgpy/data/kinetics/family.py
rmgpy/electron_placement.py
rmgpy/kinetics/arrhenius.pyx
rmgpy/reaction.py
rmgrc
test/rmgpy/data/kinetics/electronPropagationTest.py
test/rmgpy/electronPlacementTest.py
```

**[M]** Including `rmgpy/reaction.py` — the file the identity repair changes — and
`rmgpy/electron_placement.py`. If the ten had been outside the union's base, this overlap is
exactly where a merge would have conflicted.

### 1.3 But they were never outside it

```
$ git rev-list --count i123b-reaudit..plasma
0
$ git merge-base i123b-reaudit plasma
a61dc13039ff930abec5abd179675dd07ff662a6      # == plasma's own tip
```

**[M]** Local `plasma` is a strict **ancestor** of the second pass's integration branch. All ten
commits were already inside the union when the second pass built it, and are inside this pass's
union now. There was never a merge of the ten to conflict, so the scenario the brief warns
about — *"the previous passes' conflict-free readings were measured on a base that does not exist
any more"* — **cannot have happened**. The error was in the *arithmetic of the reported deltas*,
not in what was built.

### 1.4 Does the correction change either previous verdict? **No. Neither.**

Stated per pass, with the base each one used:

| pass | merge inventory taken against | unit-suite baseline taken against | verdict | changed by the correction? |
|---|---|---|---|---|
| 1 — `docs/i123-integration-audit.md` | `origin/plasma` (**stale**) | `origin/plasma`, **2596 passed / 0 failed** | NOT READY (Chemkin writer refused every plasma rate law) | **No** |
| 2 — `docs/i123b-reaudit.md` | local `plasma` | local `plasma`, **2623 passed / 0 failed** | NOT READY (sink dropped as duplicate; no Chemkin round trip) | **No** |

**[M]** Pass 1 did baseline against the stale ref — its own §2 records the divergence and its
§"Verification" table quotes `| origin/plasma baseline | 2596 | 0 | 49 | 68 |`. **Pass 2 did not.**
Its §1 graph is `git log --oneline --graph plasma..HEAD` and its §5 baseline is *"shared branch
(`plasma` + database `plasma`) … 2623 passed"* — the local ref, and the same number this pass
re-measures below. So the brief's framing (*"both prior audits measured against the stale ref"*)
is **accurate for pass 1 and inaccurate for pass 2**. What pass 2 left open was not which base it
used but which base was *authoritative*; §11 of its report lists that as an open question. The
ruling settles the question; it moves no number.

Neither verdict moves because neither rested on a delta. Pass 1's blocker was a writer that
refused a rate law — a property of the union's code, not of a count against a base. Pass 2's two
blockers were a dropped reaction and a refused read — likewise.

**[M]** The RMG-database repository has no such divergence:

```
$ git -C ../RMG-database-i123c-audit rev-list --count refs/remotes/origin/plasma..plasma ; \
  git -C ../RMG-database-i123c-audit rev-list --count plasma..refs/remotes/origin/plasma
0
0
```

Both are `fb3c13c60`. Every database number below is against that.

---

## 2. The union, and how it was built

Both fix branches descend from `93b83ac9d`, the second pass's own tip, which is also this branch's
starting point. **Both merges were conflict-free** [M]:

```
$ git merge --no-ff i134-duplicate-drops-sink   # 22 files changed, 3330 insertions(+), 5 deletions(-)
$ git merge --no-ff i135-tdep-roundtrip         # 15 files changed,  948 insertions(+), 21 deletions(-)
```

### 2.1 RMG-Py merge inventory, from git

```
$ git log --oneline --graph plasma..HEAD
*   1b7bafc45 Merge the Chemkin TDEP reader repair (I-135) into the third-pass audit union
|\
| * 948a783fd Record the I-135 evidence: before, after, red check, negative control
| * 36a43f54c Fail the run when an end-of-run export step produces nothing
| * 89c896a19 Teach the Chemkin reader the TDEP auxiliary keyword
* |   7758c27d2 Merge the reaction-identity repair (I-134) into the third-pass audit union
|\ \
| |/
|/|
| * 9b8cc38e7 Record the I-134 investigation, and the landmine the repair does not reach
| * b8204f68b Compare electron placement, not the net count, when deciding reaction identity
|/
* 93b83ac9d Report the I-123 re-audit: NOT READY, two blocking failures
* 9874fc706 Point the re-audit branch's rmgrc at the re-audit database worktree
* 1b79bee5b Let the export path read the electron-placement declaration the reactor reads
* 3559bf431 Correct the live-work exclusion section: both branches moved during the audit
* 9fc5258d9 Report the I-123 integration audit: NOT READY, one blocking failure
* ae52de38a Point the integration branch's rmgrc at the I-123 integration database
*   f087cdf04 Merge i102-quarantine: hard-fail when quarantined database data reaches a reaction model
|\
| * 99602b099 docs: evidence for the Cation_R_Recombination Marcus quarantine
| * 541e6498f kinetics: hard-fail when quarantined database data reaches a reaction model
*   5e1686bcb Merge i112-marcus-work-terms: let the Marcus work terms reach the barrier they belong in
|\
| * ab3d84072 Let the Marcus work terms reach the barrier they belong in
*   9f0eae40b Merge i115-preflight-deck: exclude the battery-SEI family from the plasma preflight decks
|\
| * ceb9fd5f1 Exclude the battery-SEI family from the plasma preflight decks
*   c1f23e93e Merge i119-rr-registry: the charge-network chain, seven tickets in one branch
|\
| * 92e2d4234 Let a kinetics library own radiative recombination
| * 18abb2789 Let a kinetics library own electron-impact ionisation
| * dde602778 Let a family declare electrons on both sides of a reaction
| * 98b0e70b4 Re-measure the electron-representation disagreement instead of arbitrating it
| * 0af07c2a2 Pin the electron-representation matrix as a measurement rather than a claim
| * 0a3b0ff3d Let the two plasma rate laws declare the electron they consume
| * 6cb64c70d Give the electron-count fixtures chemistry that actually balances
| * 04bbb2881 Propagate a library entry's declared electron count to its reaction
| * aee9e510e Make Reaction.is_balanced enforce the charge balance it computes
* 582fd6ab9 Merge i110-make-guard: refuse a bare `make` that installs into the shared conda environment
* fa042c08d Document the shared-environment guard in CLAUDE.md
* 5094dbea3 Refuse to install into the shared environment from a bare `make`
```

**Verified by inspection, not assumed** — every merge's second parent equals the branch tip it
claims [M]:

```
1b7bafc45  p2=948a783fd == i135-tdep-roundtrip
7758c27d2  p2=9b8cc38e7 == i134-duplicate-drops-sink
f087cdf04  p2=99602b099 == i102-quarantine
5e1686bcb  p2=ab3d84072 == i112-marcus-work-terms
9f0eae40b  p2=ceb9fd5f1 == i115-preflight-deck
c1f23e93e  p2=92e2d4234 == i119-rr-registry
582fd6ab9  p2=fa042c08d == i110-make-guard
```

`i126-chemkin-electrons` (`1b79bee5b`) contributes as a direct commit on the line rather than a
merge; it is present as its own tip.

### 2.2 RMG-database merge inventory, from git

```
$ git -C ../RMG-database-i123c-audit log --oneline --graph plasma..HEAD
* 2d7123b08 Give argon a sourced cation thermochemistry, and audit the gap it came from
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

All five merges' second parents match their branch tips [M]. **Neither fix branch has a database
half** — no `i134-*` or `i135-*` branch exists in RMG-database — so `i123c-audit-db` is content-
identical to `i123b-reaudit-db` at `2d7123b08`, and that is a measurement, not an omission [M].

---

## 3. The pinned-configuration hazard — the branch is **NOT merge-clean**

The tracked `rmgrc` [R] is the hazard. `fa3b58f05` (one of the ten, §1.1) created it pointing at
`../RMG-database-plasma/input` — the plasma database, which carries **neither** new kinetics
library. A suite run from a repo root with that line loads a database in which
`PlasmaElectronImpactIonization` and `PlasmaRadiativeRecombination` simply do not exist, and every
plasma assertion then passes vacuously or fails for the wrong reason. Each audit pass has had to
repoint it, and each repoint is a commit that must come out.

**There are now three such commits, not two** [M]:

| commit | pass | repointed to |
|---|---|---|
| `ae52de38a` | first | `../RMG-database-i123-integration/input` |
| `9874fc706` | second | `../RMG-database-i123b-reaudit/input` |
| `34de41c76` (this pass) | third | `../RMG-database-i123c-audit/input` |

This pass had to add a third because the second pass's line named
`../RMG-database-i123b-reaudit` — a worktree **this pass does not own**, which another branch is
free to move underneath a run in progress. Reading it would not have been wrong today; depending
on it would have been.

**Statement, plainly: the branch as it stands is NOT merge-clean.** Before any merge to `plasma`,
all three commits must be reverted so that `rmgrc` reads

```
database.directory = ../RMG-database-plasma/input
```

Nothing else on either branch needs a revert.

**Which database actually loaded, printed rather than trusted** [M] — every probe in this report
opens with the block quoted at the head of this document, and each library load logs its own
absolute path:

```
INFO:root:Loading kinetics library PlasmaElectronImpactIonization from
  /home/alon/Code/RMG-database-i123c-audit/input/kinetics/libraries/PlasmaElectronImpactIonization/reactions.py
INFO:root:Loading kinetics library PlasmaRadiativeRecombination from
  /home/alon/Code/RMG-database-i123c-audit/input/kinetics/libraries/PlasmaRadiativeRecombination/reactions.py
```

---

## 4. Build, stated before any test result

**[M]** `make build` on the union, from a worktree with **zero** `.so` files at the start:

```
$ make build
...
### make build exit: 0 ###
$ find . -name '*.so' | wc -l
104
```

The compiled modules are the ones being imported, and they carry the repair [M]:

```
compiled module = /home/alon/Code/RMG-Py-i123c-audit/rmgpy/molecule/molecule.cpython-39-x86_64-linux-gnu.so
reaction module = /home/alon/Code/RMG-Py-i123c-audit/rmgpy/reaction.cpython-39-x86_64-linux-gnu.so
rmgpy.electron_balance.get_electron_placement_counts -> <function ... at 0x...>
```

A second, incremental `make build` also exits 0.

### 4.1 The bare-`make` guard: refused on the union, **absent from the baseline**

**[M]** On the union:

```
$ make
Refusing to modify the shared RMG environment.
Use `make build` for an in-place worktree build.
Editable installation requires an explicit maintenance procedure:
    make unsafe-install-shared-env CONFIRM_SHARED_ENV_MUTATION=yes
...
make: *** [Makefile:72: guard] Error 1
```

The guard works.

**The guard is not on `plasma`.** `i110-make-guard` is one of the merges *under audit* — `582fd6ab9`
appears in `plasma..HEAD` (§2.1), not below it. By inspection [M]:

```
$ git show plasma:Makefile | grep -c 'guard\|CONFIRM_SHARED_ENV_MUTATION'
0
$ git show plasma:Makefile | grep -n '^all:'
17:all: check $(INSTALL_SENTINEL) build check
```

So on the authoritative shared branch a bare `make` still runs the editable install. **Reported as
an incident, because it is one:** I ran `make` once in the baseline worktree before establishing
this, expecting the refusal. It died at the first dependency check (SIGPIPE, exit 141) before
`pip` was reached, and the shared environment is verifiably unchanged — `PYTHONPATH= python -c
"import rmgpy"` still resolves to `/home/alon/Code/RMG-Py`, and `site-packages` carries no RMG
editable record [M]. I did not repeat the command; the Makefile inspection above is what
establishes the absence. **This is an argument for the merge, not against it.**

---

## 5. The unit suites, same command, run serially

`python -m pytest -m "not functional and not database"` — no `-n`, one run at a time, because
concurrent runs across worktrees share a fixed scratch path outside every worktree and produce
false failures.

| tree | base | failed | passed | skipped | deselected | wall |
|---|---|---|---|---|---|---|
| **union** `i123c-audit` @ `1b7bafc45` | — | **2** | **3020** | 50 | 143 | 540.78 s |
| **authoritative baseline** local `plasma` @ `a61dc1303` | — | **0** | **2623** | 50 | 83 | 503.98 s |

**[M]** Union:

```
= 2 failed, 3020 passed, 50 skipped, 143 deselected, 54 warnings in 540.78s (0:09:00) =
```

**[M]** Baseline (`/home/alon/Code/RMG-Py-i123c-step` detached at `a61dc1303`, `rmgrc` pointed at
`/home/alon/Code/RMG-database-i123c-baseline/input`, a worktree of the database at its `plasma`
tip `fb3c13c60`):

```
= 2623 passed, 50 skipped, 83 deselected, 54 warnings in 503.98s (0:08:23) =
```

The delta is **+397 passed, +2 failed, +60 deselected** against the right base. The second pass
reported 2623 on this same base, so that figure reproduces exactly.

### 5.1 The RMG-database suites

`python -m pytest test` from each database repository root, paired with the matching RMG-Py tree.

| tree | failed | passed | errors |
|---|---|---|---|
| union `i123c-audit-db` @ `2d7123b08` + union RMG-Py | **4** | **173** | 0 |
| baseline local `plasma` @ `fb3c13c60` + baseline RMG-Py | **4** | **41** | 0 |

**[M]** Union: `4 failed, 173 passed, 284 warnings in 42.85s`. Baseline:
`4 failed, 41 passed in 3.15s`. **The same four failures on both**, all in
`test/test_plasma_electron_attachment.py` and all one assertion —
`assert <KineticsDepository "Plasma_Electron_Attachment/training"> == 'rate rules'`, i.e. the O₂
attachment rate is resolving through the training depository rather than through a rate rule.
Pre-existing, identical on both bases, untouched by this union, and outside its scope. Pass 1
reported the same 4/41 baseline; it reproduces exactly.

**One finding on the way in.** Run as-is, the *baseline* database suite does not fail — it
**errors at fixture setup, 45 times**, because it resolves `database.directory` to
`/home/alon/Code/RMG-database/input`, the polymer-branch database [M]:

```
AssertionError: RMG resolved database.directory to '/home/alon/Code/RMG-database/input',
not this worktree ('/home/alon/Code/RMG-database-i123c-baseline/input').
```

The baseline has no `test/conftest.py`; the union does, because `bdbbf2864` *"test: pin
database.directory for the suites in this repository"* is one of the merges under audit (§2.2).
The baseline number above was taken with an equivalent `rmgrc` written at the database repository
root, which is the same workaround pass 1 had to use. **That pin is a union improvement, and it is
the difference between a suite that measures its own tree and one that silently measures another.**

---

## 6. The two union failures, named and attributed

Both are **pre-existing on the union's own base and neither is union-only** [M]:

```
FAILED test/rmgpy/electronPlacementTest.py::TestElectronPlacementResolver::test_declaration_registry_is_explicit_and_closed
  AssertionError: assert {'Plasma_Electron_Attachment': (1, 0), 'Cation_R_Recombination': (1, 0),
                          'PlasmaRadiativeRecombination': (1, 0), 'PlasmaElectronImpactIonization': (1, 2)}
                      == {'Plasma_Electron_Attachment': (1, 0), 'Cation_R_Recombination': (1, 0),
                          'PlasmaElectronImpactIonization': (1, 2)}

FAILED test/rmgpy/preflightDeckFamilyExclusionTest.py::PlasmaDeckFamilyExclusionTest::test_plasma_deck_does_not_declare_family[docs/i102-quarantine/input.py]
  AssertionError: plasma deck docs/i102-quarantine/input.py declares Cation_R_Recombination, a
  lithium-ion-battery SEI family excluded from every plasma configuration.
```

1. **The registry closure assertion.** `i119-rr-registry` added `PlasmaRadiativeRecombination` to
   `FAMILY_ELECTRON_PLACEMENT` [R] `rmgpy/electron_placement.py` and did not update the test that
   asserts the registry's exact contents. The test asserts three entries against the shipped four.
   **The code is right and the test is stale.** Attributed to `i119-rr-registry` alone — see §7,
   where it is red on that branch by itself.
2. **The quarantine deck declaring the quarantined family.** `docs/i102-quarantine/input.py` exists
   *to* declare `Cation_R_Recombination` — that is what it reproduces — and `i115-preflight-deck`'s
   sweep test globs every plasma deck in the tree, including that one. Two branches whose
   individual intents are both correct; the interaction is a test-discovery scope question.
   Attributed to the pair **`i102-quarantine` × `i115-preflight-deck`**.

**Union-only failures: there are none.** Every failure the union shows is reproduced on a
contributing branch or on a contributing pair; §7 is the measurement. Neither is a wrong number
and neither is mine to fix here — the brief forbids fixing a failing test by changing the test,
and #1 is precisely a test that needs changing by whoever owns `i119`.

---

## 7. Each contributing branch's own state, measured this pass

Measured, not inherited, in a single stepping worktree
(`/home/alon/Code/RMG-Py-i123c-step`) checked out detached at each tip in turn, `check-pydas`'d,
built in place and run with the same command — **strictly serially**.

| branch | tip | database it was paired with | build | failed | passed | skipped | deselected |
|---|---|---|---|---|---|---|---|
| `i110-make-guard` | `fa042c08d` | baseline `fb3c13c60` | exit 0 | **0** | 2635 | 50 | 68 |
| `i119-rr-registry` | `92e2d4234` | baseline `fb3c13c60` | exit 0 | **1** | 2684 | 50 | 118 |
| `i115-preflight-deck` | `ceb9fd5f1` | baseline `fb3c13c60` | exit 0 | **0** | 2630 | 50 | 87 |
| `i112-marcus-work-terms` | `ab3d84072` | baseline `fb3c13c60` | exit 0 | **0** | 2667 | 50 | 83 |
| `i102-quarantine` | `99602b099` | baseline `fb3c13c60` | exit 0 | **0** | 2641 | 49 | 90 |
| `i126-chemkin-electrons` | `1b79bee5b` | first-pass `c7bd96292` | exit 0 | **2** | 2835 | 49 | 131 |
| `i134-duplicate-drops-sink` | `9b8cc38e7` | union `2d7123b08` | exit 0 | **2** | 3002 | 50 | 143 |
| `i135-tdep-roundtrip` | `948a783fd` | union `2d7123b08` | exit 0 | **2** | 2853 | 49 | 131 |
| *(union, for comparison)* | `1b7bafc45` | union `2d7123b08` | exit 0 | **2** | 3020 | 50 | 143 |

Each branch was paired with the database its **own** tracked `rmgrc` names, remapped onto a
worktree this pass owns so that no other worker's tree is depended on: `../RMG-database-plasma`
→ `RMG-database-i123c-baseline` (`fb3c13c60`), `../RMG-database-i123-integration` →
`RMG-database-i123c-step1` (`c7bd96292`), `../RMG-database-i123b-reaudit` →
`RMG-database-i123c-audit` (`2d7123b08`). `i110-make-guard` carries **no** tracked `rmgrc` at all
— re-confirming the second pass's finding — and was pinned to the baseline database explicitly.

**Green/red, measured this pass and inherited from nobody:**

- **Red on its own branch: `i119-rr-registry`, one failure** — the registry-closure assertion
  (§6.1). This is the branch the second pass found red, and it is still red. **[M]**
  ```
  = 1 failed, 2684 passed, 50 skipped, 118 deselected, 54 warnings in 617.51s (0:10:17) =
  FAILED test/rmgpy/electronPlacementTest.py::TestElectronPlacementResolver::test_declaration_registry_is_explicit_and_closed
  ```
- **Green on their own branches: `i110-make-guard`, `i115-preflight-deck`,
  `i112-marcus-work-terms`, `i102-quarantine`** — zero failures each. In particular
  **`i102-quarantine` and `i115-preflight-deck` are each green alone**, which is what makes §6.2
  an interaction rather than a defect in either.
- **Red on their own branches with the union's exact pair of failures: `i126-chemkin-electrons`,
  `i134-duplicate-drops-sink`, `i135-tdep-roundtrip`.** All three sit above the merges of both
  `i119` and the `i102`/`i115` pair, so they inherit both.

**Union-only failures: there are none.** Every failure the union shows is already present on a
contributing branch — the first from `i119-rr-registry` alone, the second from `i126` downward,
i.e. from the first point at which `i102-quarantine` and `i115-preflight-deck` are both in the
tree. Nothing in the merged suite fails that passes on every contributor.

`python utilities.py check-pydas` was run explicitly before each `make build`, because the
`Makefile` on the older branches has no `check-pydas` prerequisite on `build` and a fresh worktree
has no `rmgpy/solver/settings.pxi`. That deviation is the same one the second pass recorded; it is
a harness step, not a code change, and it is why every row above reads *build exit 0*.


---

## 8. THE BLOCKING FAILURE — a re-exported plasma mechanism is wrong by 10⁶

### 8.1 What happens

`scripts/mergeModels.py`, unmodified, on two Chemkin models each carrying one half of the lithium
charge network — the files produced by the deck run in §10, split one reaction each [M]:

```
$ python scripts/mergeModels.py --model1 .../merge-half0/chem.inp .../species_dictionary.txt \
                                --model2 .../merge-half1/chem.inp .../species_dictionary.txt
Added 1 out of 1 (100.0%) unique reactions from model #1.
Added 1 out of 1 (100.0%) unique reactions from model #2.
The merged model has 7 species and 2 reactions
### mergeModels exit: 0 ###
```

Both channels survive the merge — that part is the I-134 repair working. The numbers do not:

```
                                                    A          n        Ea
as generated  Li(2)+e-(1)=>[Lip](3)+e-(1)+e-(1)     1.490e+18  -0.267   162.150
after merge   Li(2)+e-(1)=>[Lip](3)+e-(1)+e-(1)     1.490e+12  -0.267   162.150
as generated  [Lip](3)+e-(1)=>Li(2)                 1.734e+14  -0.801     0.010
after merge   [Lip](3)+e-(1)=>Li(2)                 1.734e+08  -0.801     0.010
```

Both files carry the identical header `REACTIONS    KCAL/MOLE   MOLES` [M]. Exit status 0, no
warning, no log line. A mechanism that is 10⁶ slow in both its ionisation source and its
recombination sink.

### 8.2 Isolated to a read-then-write, without `mergeModels` in the way

`load_chemkin_file` then `save_chemkin_file`, nothing else [M]
(`evidence/probe-rewrite-units.log`):

```
read back:
  Li(2) + e-(1) => [Lip](3) + e-(1) + e-(1)
    kinetics : TwoTemperaturePlasma
    A        : 1.49e+18 cm^3/(mol*s)   (value_si = 1.490000e+12)
    k(11604.5 K, Te=11604.5 K) = 1.081914e+08
as written the SECOND time:
  Li(2)+e-(1)=>[Lip](3)+e-(1)+e-(1)                   1.490e+12 -0.267    162.150
```

**The reader is correct.** The reloaded object carries `A = 1.49e+18 cm^3/(mol*s)`, whose SI value
is `1.49e+12`, and reproduces the generated rate to every printed digit. **The writer is what is
wrong.**

And it compounds: reading the second file gives a rate 10⁶ below the first [M]:

```
Li(2) + e-(1) => [Lip](3) + e-(1) + e-(1)
    k after one round trip  = 1.081914e+08
    k after two round trips = 1.081914e+02   ratio 1.000000e-06
[Lip](3) + e-(1) => Li(2)
    k after one round trip  = 9.618226e+04
    k after two round trips = 9.618226e-02   ratio 1.000000e-06
```

### 8.3 The cause, one line

**[R]** `rmgpy/chemkin.pyx`, `_plasma_arrhenius_for_chemkin`:

```python
if isinstance(kinetics, _kinetics.TwoTemperaturePlasma):
    arrhenius = _kinetics.Arrhenius(
        A=(kinetics.A.value_si, kinetics.A.units),      # <-- SI VALUE, DECLARED UNITS
        n=kinetics.n.value_si,
        Ea=(kinetics.Ea_g.value_si, 'J/mol'),           # value_si with SI units: correct
        T0=(kinetics.T0.value_si, 'K'),                 # value_si with SI units: correct
    )
```

An SI *value* paired with whatever units string the object happens to *declare*. Driven directly on
two objects that are the same rate law written two ways [M]:

```
SI  : A = 1.49e+12 m^3/(mol*s)    value_si = 1.490000e+12
cgs : A = 1.49e+18 cm^3/(mol*s)   value_si = 1.490000e+12      <- the same rate law
reduced for Chemkin:
  from the SI-declared object  : A = 1.49e+12 m^3/(mol*s)
  from the cgs-declared object : A = 1.49e+12 cm^3/(mol*s)     <- 10^-6 of the right number
  ratio of the two reductions  : 1.000000e-06
```

Ten to the minus six for a second-order rate; ten to the minus twelve for a third-order one.

### 8.4 Why it is only now reachable — and where it belongs

The writer line is **older than the union**. It came in with `a9c81e6de` "Export plasma reactions
instead of dropping them, and carry the electron", which is an ancestor of local `plasma` [M]:

```
$ git log --oneline -S "A=(kinetics.A.value_si, kinetics.A.units)" -- rmgpy/chemkin.pyx
a9c81e6de Export plasma reactions instead of dropping them, and carry the electron
$ git merge-base --is-ancestor a9c81e6de plasma && echo IN-PLASMA
IN-PLASMA
```

It was harmless there because **every `TwoTemperaturePlasma` RMG builds for itself declares SI
units**. All three producers pass `A=(arr.A.value_si, "m^3/(mol*s)")` [R]
`rmgpy/kinetics/arrhenius.pyx:802, 1428, 2037`, and no database entry anywhere uses the class
[D] (`grep -rl TwoTemperaturePlasma` over the database input tree: no hits). For an SI-declared
object the wrong pairing is accidentally right.

**The Chemkin reader added by I-135 is the one producer of a `TwoTemperaturePlasma` with non-SI
units** [R] `rmgpy/chemkin.pyx::_two_temperature_plasma_from_arrhenius`, which rebuilds the rate
law from the file's own `cm^3/(mol*s)` numbers: `A=(arrhenius.A.value, arrhenius.A.units)`.

On the authoritative baseline the whole path is unreachable, because the read refuses [M]:

```
$ PYTHONPATH=/home/alon/Code/RMG-Py-i123c-baseline python -c "load_chemkin_file(...)"
rmgpy: /home/alon/Code/RMG-Py-i123c-baseline/rmgpy/__init__.py
BASELINE REFUSED: ChemkinError 'e-(1)' doesn't look like a collision efficiency for species TDEP in line 'TDEP/E-(1)/'
```

**So: old defect, new reachability.** Not a defect introduced by the merge, and not one the merge
may ship either — a loud refusal has been replaced by a silent wrong number, which is the exact
trade this campaign has repeatedly refused. It is not union-only in the narrow sense: it is
reachable on `i135-tdep-roundtrip` by itself, since that branch carries both halves — **measured,
not inferred**, by running the same probe in the stepping worktree checked out detached at
`948a783fd` [M]:

```
rmgpy: /home/alon/Code/RMG-Py-i123c-step/rmgpy/__init__.py     # HEAD 948a783fd i135-tdep-roundtrip
A column, first write : [1.49e+18, 1.734e+14]
A column, second write: [1.49e+12, 1.734e+08]
    1.490000e+18 -> 1.490000e+12   ratio 1.000000e-06
    1.734000e+14 -> 1.734000e+08   ratio 1.000000e-06
```

### 8.5 Why this blocks

- It is a **laundered quantity**: correct in, wrong out, exit 0, no warning, in a file whose header
  states units the numbers are not in.
- It is reachable by a **shipped script** (`scripts/mergeModels.py`) on the campaign's own deck
  output, with no plasma-specific invocation.
- The read-back capability the union adds is what makes it reachable, so it arrives *with* the
  merge rather than being inherited by it.
- The fix is small and local — pair `value_si` with SI units, as the same constructor already does
  for `Ea` and `T0` two lines down — but it is **not mine to make**: this pass audits, and the
  brief forbids fixing.

---

## 9. The three previously-unexercised consumers of the repaired identity predicate

`docs/i134-duplicate-electrons.md` §9 marks all three UNKNOWN. All three are now measured
(`evidence/probe-identity-consumers.log`).

### 9.1 Model merging — **exercised, correct, with both controls**

`ReactionModel.merge` [R] `rmgpy/rmg/model.py:129`, consuming
`is_isomorphic(either_direction=True)`. Two models, one channel each, the real library reactions
[M]:

```
ionisation    : [Li] => [Lip]  electrons=1   placement=(1, 2)  owner=PlasmaElectronImpactIonization
recombination : [Lip] => [Li]  electrons=-1  placement=(1, 0)  owner=PlasmaRadiativeRecombination

merged model: 2 species, 2 reactions
  - [Li] => [Lip]   electrons=1   kinetics=VoronovEIArrhenius
  - [Lip] => [Li]   electrons=-1  kinetics=BadnellRRArrhenius

[PASS] merge keeps BOTH channels -- 2 reaction(s) in the merged model
[PASS] the merged model carries both rate laws
[PASS] the cation is consumed somewhere in the merged model
```

Both controls hold, so "both survived" does not merely mean merge stopped deduplicating [M]:

```
[PASS] CONTROL A: merging a model with itself still deduplicates      -- 1 reaction(s)
[..]   true reverse of the ionisation: [Lip] => [Li]  electrons=-1  placement=(2, 1)
[PASS] CONTROL B: a genuine reverse pair still collapses in merge     -- 1 reaction(s)
```

**No other electron-blind step was found in merge.** The end-to-end `mergeModels.py` run in §8.1
also keeps both channels — and is where the separate rate defect surfaced.

### 9.2 The pressure-dependence path — **latent, and measured rather than assumed**

**[M]** No shipped plasma deck turns pressure dependence on. Five plasma decks in the tree
(`docs/i102-quarantine/input.py`, `docs/i123-integration/input.py`, `docs/m7-preflight/input.py`,
`input_placement.py`, `input_secondary.py`); none contains a `pressureDependence(` block.

**[M]** Nothing in `rmgpy/rmg/pdep.py` ever gives a `PDepReaction` a non-zero electron count — the
three mentions of `electrons` in the whole file are the constructor default `electrons=0`, the
pass-through, and `__reduce__`:

```
[PASS] a PDepReaction is born with electrons = 0 -- 0
[PASS] and therefore with placement counts (0, 0) -- (0, 0)
```

**[M]** The comparison `pdep.py:927/939` actually makes, driven directly:

```
LibraryReaction ionisation (electrons=+1)     vs neutral PDepReaction -> is_isomorphic(either_direction=True) = False
LibraryReaction recombination (electrons=-1)  vs neutral PDepReaction -> is_isomorphic(either_direction=True) = False
```

**[I]** Basis: a `PDepReaction` is necessarily neutral, so its placement counts are `(0, 0)`; a
plasma library reaction's are `(1, 2)` or `(1, 0)`. They cannot match, before or after the repair —
the old net-scalar test also separated `±1` from `0`. **So the repair changes no verdict here, and
the path is genuinely latent.** Stated as such rather than claimed as coverage.

### 9.3 The products-only generation shortcut — **reachable; the doubted negation is not**

`Reaction.is_isomorphic(check_template_rxn_products=True)` [R] `rmgpy/reaction.py:707-727`, which
still compares the net `electrons` scalar and negates it on `is_forward`. The question the brief
poses is scoping.

**It is reachable.** Its only production callers are inside `find_degenerate_reactions`
[R] `rmgpy/data/kinetics/common.py:356,359` and its family-side siblings
[R] `rmgpy/data/kinetics/family.py:1950,1963`. Driving the real route RMG itself takes —
`KineticsDatabase.generate_reactions_from_families` with `Plasma_Electron_Attachment` on triplet
O₂, out of the union database, no lithium in it [M]:

```
1 reaction(s) generated by Plasma_Electron_Attachment from triplet O2, through
generate_reactions_from_families -- the route RMG itself takes
  - [O][O] => [O][O-]   electrons=-1  is_forward=True  degeneracy=1.0  type=TemplateReaction

find_degenerate_reactions calls recorded during that generation:
  call 0: offered=2 returned=1 charged_electrons=[-1, -1] is_forward=[True, True] comparisons_possible=True

[PASS] the degeneracy path -- the shortcut's only production caller -- IS reached during plasma family generation
[PASS] and it is reached carrying CHARGED reactions
[PASS] and with more than one reaction offered, so a comparison actually happens
```

Two charged reactions offered, one returned: a real deduplication, on the shortcut, with charged
reactions.

**The doubted negation is not reachable on that path.** Every reaction entering the shortcut here
carries `is_forward=True`, so `electrons1 = self.electrons` and the `-self.electrons` branch is
never taken [M]:

```
[PASS] every reaction reaching the shortcut on this path is is_forward=True, so the
       doubted negation is never taken here -- {True}
```

And the shortcut does separate a charged reaction from its charge-reversed twin [M]:

```
same electrons     -> True
opposite electrons -> False
```

**Answer: reachable, and correct on what it reaches. The negation whose correctness I-134's author
doubted is not exercised by any plasma generation path measured here.** Per the brief this is where
it stops — nothing was changed. **UNKNOWN** remains whether some other caller can present an
`is_forward=False` reaction to it; `family.py:1950,1963` were not driven.

---

## 10. The three known defects, re-confirmed on this union

`evidence/probe-known-defects.log`. Confirmed as described, not re-litigated, not fixed.

### 10.1 The RMS YAML writer still refuses these rate laws — and still **loudly**

**[M]** The property the decision rests on is that it raises at save rather than emitting a wrong
number, and it does:

```
rmgpy/yaml_rms.py:252  raise ValueError("Object of type {} does not have a defined conversion to "
  preceding line: else:
[PASS] the dispatch's terminal branch is a raise, not a fall-through default

PlasmaElectronImpactIonization -> VoronovEIArrhenius  LOUD: ValueError
    Object of type <class 'rmgpy.kinetics.arrhenius.VoronovEIArrhenius'> does not have a
    defined conversion to ReactionMechanismSimulator format
PlasmaRadiativeRecombination   -> BadnellRRArrhenius  LOUD: ValueError
[PASS] both rate laws are REFUSED loudly by the RMS writer
```

**[D][M]** The deck-level workaround is still in place: `docs/i123-integration/input.py` sets
`generateRMSYAML=False`, and the deck run in §11 wrote no `rms/` directory while the non-plasma
negative control did. Kinetics-coverage gap, not an electron-representation one. **Unchanged, not
blocking.**

### 10.2 One library carrying both channels still cannot be loaded

**[M]** The shipped entries of the two libraries, merged into one library and offered to
`check_for_duplicates`:

```
entry 1: [Li] => [Lip]   electrons=1   family=None  placement=(0, 1)
entry 2: [Lip] => [Li]   electrons=-1  family=None  placement=(1, 0)
check_for_duplicates REFUSED -> DatabaseError
    Unexpected duplicate reaction [Li] => [Lip] in kinetics library PlasmaBothChannels.
    Reaction index 2 matches index 1.
```

`family=None` on both is the cause and it is upstream of the repair: `load_entry` builds a plain
`Reaction` with no owner, so no declaration is visible and the net rule makes them exact mirrors.
The owner is attached later, in `get_library_reactions`.

**[M]** The pin is still expected-failing rather than silently passing: the strict
`@pytest.mark.xfail` on `test_one_library_carrying_both_channels_can_be_loaded`
(`test/rmgpy/i134DuplicateElectronsTest.py`) is in the union suite's **50 skipped / xfailed**
column and does not appear among the 2 failures — a strict xfail that started passing would show
as `XPASS`-failure. **[M]** The workaround holds: the shipped deck still declares
`PlasmaElectronImpactIonization` **and** `PlasmaRadiativeRecombination` as two separate
`reactionLibraries` entries, and each library carries exactly one entry. **Unchanged, not
blocking.**

### 10.3 The lithium cation's enthalpy — **measured on this union, not inherited**

**[D][M]** Straight out of the database this union loads:

```
lithium/cation thermo libraries present: ['LithiumPrimaryThermo', 'LithiumAdditionalThermo', 'PlasmaCationThermo']
LithiumPrimaryThermo :: entry '[Lip]' index=65
    H298 = 532.936 kJ/mol
    E0   = 526.738 kJ/mol
```

**[I]** Reference: `dfH(Li⁺,g) = dfH(Li,g) + IE(Li) = 159.3 + 520.2 = 679.5 kJ/mol` — standard
atomisation enthalpy and first ionisation energy of lithium, quoted rather than re-derived here.

```
LithiumPrimaryThermo :: [Lip]  H298 = 532.936 kJ/mol   gap vs reference = +146.6 kJ/mol (+1.52 eV)
```

**The defect is still present on this union, at the same ~1.5 eV.** The database branch that
diagnosed it (`i129-lithium-cation-enthalpy`, tip `75e0cb6f0`) is **not** merged into
`i123c-audit-db` (§2.2), so the union carries the old value. The deck confirms it uses this entry:
`chem_annotated.inp` writes `[Lip](3)` with `a6 = 6.33520000E+04`, i.e. `H/R = 63352 K`,
`H = 526.7 kJ/mol` — the `E0` above [M].

Non-blocking for the same reason as before: the two reactions in the deck are irreversible, so no
`Keq` is formed from this enthalpy and no rate in the generated model depends on it. **But it is
not fixed, and this union is not the tree where it is fixed.**

---

## 11. The lithium charge network, end to end, stage by stage

`evidence/probe-lithium-network.log`, `evidence/probe-tdep-readback.log`,
`evidence/probe-readback-identity.log`, and the deck run.

| # | stage | result |
|---|---|---|
| 1 | **Load** both kinetics libraries out of the union database | PASS |
| 2 | **Balance** each canonical reaction | PASS |
| 3 | **Placement** resolution from the declaration | PASS |
| 4 | **Reactor acceptance** — one `PlasmaReactor` holding source and sink | PASS |
| 5 | **Presence in the generated model** — not merely loadable | PASS |
| 6 | **Write-out** | PASS |
| 7 | **Read-back** | PASS, with the representation property in §11.3 |
| 8 | **Re-export after read-back** | **FAIL — §8** |

### 11.1 Stages 1–4, from the loader, the resolver and the reactor

**[M]**, all twenty-two checks passing:

```
STAGE 1 -- LOAD both kinetics libraries out of the integration database
  [PASS] ionisation library loaded / recombination library loaded
  [PASS] ionisation declares a net electron GAIN (+1) / recombination a net LOSS (-1)
  [PASS] the neutral feed and the cation are the same two species in both channels
STAGE 2 -- BALANCE
  [PASS] ionisation is_balanced / recombination is_balanced
STAGE 3 -- PLACEMENT
  [PASS] both libraries are declared owners; ionisation (1, 2), recombination (1, 0)
  [PASS] ionisation resolves to Li + e- => Li+ + 2 e-
  [PASS] recombination resolves to Li+ + e- => Li + hv
  [PASS] both views carry no residual electron metadata
  [PASS] resolution did not mutate the canonical reactions
STAGE 4 -- REACTOR
  [PASS] the reactor accepted both reactions (2) and located the electron
  [PASS] source evaluated at the Voronov rate / sink evaluated at the Badnell rate
  [PASS] the cation has a non-zero loss channel
```

### 11.2 Stage 5 and 6 — presence in the model, and the file

**[M]** `python rmg.py <deck>/input.py`, the shipped deck through the real input-file path:

```
Adding reaction library PlasmaElectronImpactIonization to model edge...
Adding reaction library PlasmaRadiativeRecombination to model edge...
MODEL GENERATION COMPLETED
The final model core has 7 species and 2 reactions
```

**Two** core reactions, and **no `This library reaction was not new:` line anywhere in `RMG.log`**
[M] — the second pass's blocker is gone. `chemkin/chem.inp`, `REACTIONS` block in full:

```
REACTIONS    KCAL/MOLE   MOLES

Li(2)+e-(1)=>[Lip](3)+e-(1)+e-(1)                   1.490e+18 -0.267    162.150
    TDEP/e-(1)/   ! VoronovEIArrhenius exported as a modified-Arrhenius fit of k(Te) over 1.16e+04-2.321e+08 K; ...

[Lip](3)+e-(1)=>Li(2)                               1.734e+14 -0.801    0.010
    TDEP/e-(1)/   ! BadnellRRArrhenius exported as a modified-Arrhenius fit of k(Te) over 10-1e+07 K; ...

END
```

`chem_annotated.inp` attributes them to `PlasmaElectronImpactIonization` and
`PlasmaRadiativeRecombination` respectively [M]. The cation appears on a reactant side: it has a
loss channel.

**The run exits 1, by design** [M]:

```
Error: EXPORT FAILED. 1 of RMG's output steps did not produce their files:
Error:   Chemkin-to-Cantera translation (cantera_from_ck/)
Error:       InputError: Error while reading reaction in chem.inp starting on line 66:
### RMG_EXIT=1 ###
```

That is I-135's exit-status change doing what it was written to do: `cantera_from_ck/` is empty
and the run says so instead of exiting 0. `cantera2/` — the direct writer, the usable Cantera
artifact — carries six files [M]. **Stated so nobody is surprised later: on this union every
plasma deck exits non-zero at the last step**, and any automation that treats a non-zero exit as
"the run failed" will need the distinction. That is a consequence of the merge, correctly
reported by it; it is not the blocking failure.

### 11.3 Stage 7 — the read-back, and what the representation means

**[M]** RMG's reader on the file as written, and the two controls:

```
A  RMG reader, file as written    : READ  -- 7 species, 2 reactions, kinetics TwoTemperaturePlasma,
                                     uses Te = True, k(11604.5 K) = 1.081914e+08 (the generated rate)
B  RMG reader, TDEP line removed  : READ  -- kinetics Arrhenius, uses Te = False
C  Cantera ck2yaml, file as written: REFUSED -- could not convert string to float: 'e-(1)'
```

Trip B is the control that isolates the `TDEP/` line as the whole of the obstacle, and it also
shows what the repair buys: without it the same numbers come back as a **gas-temperature**
Arrhenius, which is a different reaction, not a lossier one.

**The representation property, established rather than glossed** [M]:

```
canonical original : [Li] => [Lip]                                 electrons=+1  placement=(1,2)
                     reactants [('[Li]', is_electron=False)]        owner=PlasmaElectronImpactIonization
reloaded           : Li(2) + e-(1) => [Lip](3) + e-(1) + e-(1)     electrons=0   placement=(0,0)
                     reactants [('Li', False), ('e-', True)]        owner=Unclassified
```

**Q1 — does a reloaded reaction compare as the same reaction as its canonical original, under the
newly repaired predicate? NO, in either direction** [M]:

```
PlasmaElectronImpactIonization vs reloaded ionisation    same-direction=False either-direction=False
PlasmaElectronImpactIonization vs reloaded recombination same-direction=False either-direction=False
PlasmaRadiativeRecombination   vs reloaded ionisation    same-direction=False either-direction=False
PlasmaRadiativeRecombination   vs reloaded recombination same-direction=False either-direction=False
[PASS] NO reloaded reaction compares equal to its canonical original
```

**As a property, not a defect:** the round trip recovers the *chemistry* in the expanded electron
representation, not the canonical *object*. A reloaded mechanism and the model that produced it are
not interchangeable under RMG's identity predicate. Anyone who round-trips a plasma mechanism and
then compares it against the generated one — a regression diff, a model-comparison script — will
get "different" for every reaction, and that is the representation talking, not the chemistry.

**Q2 — do the two reloaded channels collapse into each other, the I-134 defect seen from the other
side of the file? No** [M]:

```
Li(2) + e-(1) => [Lip](3) + e-(1) + e-(1)  vs  [Lip](3) + e-(1) => Li(2)
    -> is_isomorphic(either_direction=True) = False
```

In the expanded form the electrons are ordinary participants, so the two channels differ in their
participant lists alone. **The expansion is self-distinguishing and needs no placement
declaration** — which is a genuinely reassuring result: the file format does not depend on the
registry the canonical form does.

**Q3 — can a reloaded plasma mechanism be fed back into a model? Two separate answers, separated by
a neutral control.**

*Charge balance — plasma-specific, and this audit's finding* [M]:

```
Li(2) + e-(1) => [Lip](3) + e-(1) + e-(1)
    reactant elements {'Li': 1, 'e': 1}   product elements {'Li': 1, 'e': 2}   electrons scalar 0
    is_balanced -> False
[Lip](3) + e-(1) => Li(2)
    reactant elements {'Li': 1, 'e': 1}   product elements {'Li': 1}   electrons scalar 0
    is_balanced -> False
```

`Reaction.is_balanced` counts the electron as a conserved *element* and then folds the scalar
`electrons` into the net charges [R] `rmgpy/reaction.py:1688-1760`. The expanded form puts the
electron in the participant lists **and** sets the scalar to zero, so the element tallies differ
and the scalar cannot compensate. `KineticsLibrary.load` raises `DatabaseError` on any entry that
fails this check [R] `rmgpy/data/kinetics/library.py:585` — **so a plasma mechanism that has been
round-tripped through Chemkin cannot be handed back to RMG as a kinetics library or a seed
mechanism.** The neutral control separates this cleanly: all 66 reactions of the non-plasma
mechanism reload with `is_balanced -> True` [M].

Not blocking, because the route RMG itself uses for restart does not go through the Chemkin file:
the deck's own `seed/seed/reactions.py` is written in the **canonical** form — `label = "Li =>
[Lip]"`, `kinetics = VoronovEIArrhenius(..., electrons=1, ...)`, no electron participants [M]. The
seed round trip is unaffected.

*Owner loss — general, and not this campaign's* [M]:

```
reloaded plasma  : owner labels ['Unclassified'] -> make_new_reaction raises
                   "Could not retrieve the family/library: Unclassified"
reloaded NEUTRAL : owner labels ['Unclassified'] -> make_new_reaction raises, identically, 5 of 5
```

Both `chem.inp` and `chem_annotated.inp` reload as `Unclassified`. `CoreEdgeReactionModel.make_new_reaction`
resolves a family/library object for every reaction [R] `rmgpy/rmg/model.py:461` via
`get_family_library_object` [R] `:2308`, which raises on an unknown label. That function is
pre-existing on `plasma` and upstream of it [M], and the neutral control fails the same way, so
**this is a general property of any Chemkin-reloaded RMG mechanism and is not an
electron-representation matter.**

### 11.4 The ecosystem Chemkin converter — the smaller claim, said small

**[M]** Cantera's `ck2yaml` (cantera 3.1.0) does **not** implement the `TDEP` auxiliary keyword and
refuses the file:

```
Error while reading reaction in chem.inp starting on line 66:
"""
Li(2)+e-(1)=>[Lip](3)+e-(1)+e-(1)                   1.490e+18 -0.267    162.150
    TDEP/e-(1)/   ! VoronovEIArrhenius exported as ...
"""
could not convert string to float: 'e-(1)'
```

**The claim this union can make is the smaller one: RMG can now read back the plasma Chemkin file
RMG wrote. The file still does not import into the Cantera ecosystem, and the union does not change
that.** The usable Cantera artifact for a plasma mechanism is `cantera2/`, written by RMG's direct
writer, which is unaffected.

---

## 12. The negative control

`examples/rmg/minimal`, run on the union and on the authoritative baseline, same input file [M]:

| | union `1b7bafc45` | baseline `a61dc1303` |
|---|---|---|
| exit status | **0** | **0** |
| core | 26 species / 66 reactions | 26 species / 66 reactions |
| edge | 146 species / 332 reactions | 146 species / 332 reactions |

**[M]** `diff -r -q` over every generated artifact, excluding `*.log`, `*.png`, `*.xls`, `*.h5` and
the copied input:

```
Files .../nc-baseline/cantera_from_ck/chem_annotated.yaml and .../nc-union/... differ
Files .../nc-baseline/cantera_from_ck/chem.yaml            and .../nc-union/... differ
```

Two files, and the entire diff of each is one line:

```
4c4
< date: Thu, 27 Aug 2026 10:07:05 +0300
---
> date: Thu, 27 Aug 2026 09:54:54 +0300
```

The conversion timestamp. Every Chemkin file, both species dictionaries, `tran.dat`, every RMS
file, the seed mechanism and every solver CSV are **byte-identical**, and no file exists on one
side only. **An ordinary non-plasma mechanism builds, runs and exports byte-identically, with the
same exit status as the authoritative shared branch.**

---

## 13. Verdict

**NOT READY.**

**The blocking failure: §8.** A plasma Chemkin mechanism that RMG reads back and writes out again
carries every rate constant 10⁶ too small, silently, under a units header the numbers are not in,
with exit status 0. Reachable through the shipped `scripts/mergeModels.py` on the campaign's own
deck output. The writer line is older than the union; the reader that makes it reachable is one of
the two repairs this pass merged, so the defect arrives *with* the merge. One line in
`rmgpy/chemkin.pyx` — `A=(kinetics.A.value_si, kinetics.A.units)` should pair `value_si` with SI
units, as the `Ea` and `T0` arguments two lines below it already do — plus a test that round-trips
a plasma mechanism twice and compares `k`.

**What is genuinely fixed, and should not be re-litigated:**

- the cation's loss channel is no longer discarded as a duplicate of its own source; both channels
  reach the generated model, in both offer orders, through the real input-file path (§11.2);
- `ReactionModel.merge` keeps both channels, with both controls holding (§9.1);
- the Chemkin file RMG writes for a plasma mechanism reads back into RMG, with its
  electron-temperature dependence intact (§11.3);
- a failed export step no longer reports success in the exit status (§11.2);
- the negative control is byte-identical (§12).

**Also required before any merge, independently of the blocker:** the three `rmgrc` pin commits
must be reverted (§3). And whoever owns `i119-rr-registry` owns the stale registry-closure
assertion (§6).

---

## 14. What this pass could not reach

- **Two stray files.** The unit suite leaves `chem-gas.yaml` and `chem_annotated.yaml`
  untracked in the repository root — `ck2yaml` conversions of the surface-chemistry test fixtures,
  written relative to the current working directory rather than into a temporary one. **[I]** Not
  union-introduced: `test/rmgpy/yaml_cantera2Test.py`, `test/rmgpy/cantera_yaml_comparer.py` and
  `test/rmgpy/rmg/mainTest.py` are byte-identical between `plasma` and the union
  (`git diff plasma..HEAD --stat` over the three: empty). I did **not** observe the baseline tree
  produce them directly, because the branch sweep's `git clean` removed untracked files between
  runs, so the attribution rests on the diff rather than on a sighting. Cosmetic; noted so the
  next reader does not mistake them for audit output.
- **Functional and database-marked RMG-Py tests.** Not run — the brief asked for the unit suite,
  and `-m "not functional and not database"` is what was measured on every tree.
- **Regression tests.** `test/regression/` is outside pytest and was not run on any tree.
- **`family.py:1950,1963`** — the two family-side callers of the products-only shortcut were
  identified statically but not driven. Whether either can present an `is_forward=False` reaction
  to it is **UNKNOWN** (§9.3).
- **Whether the 10⁶ defect also reaches the Cantera or RMS writers.** Only the Chemkin writer was
  measured. Both other writers refuse or re-derive rather than re-emit a reloaded object, so the
  shape differs — but it was not measured, and a sibling audit of the same constructor family is
  reported elsewhere in this campaign. **UNKNOWN.**
- **Performance.** The identity repair adds an O(1) comparison to RMG's hottest path. Union suite
  540.8 s over 3020 tests versus baseline 504.0 s over 2623 — 0.179 s/test versus 0.192 s/test,
  which is not a controlled benchmark on a quiet machine and should not be read as one. **No
  material slowdown observed; not measured properly.**
- **A physical validation of the lithium network.** The deck runs and both channels evaluate; that
  the resulting number densities are physically right is a separate milestone's question, and the
  known 1.5 eV cation enthalpy error (§10.3) is one reason to expect they are not.
- **Whether `origin/plasma` should be fast-forwarded.** Out of scope: the pass may not push.

---

## 15. Reproducing this

```bash
cd /home/alon/Code/RMG-Py-i123c-audit
conda activate rmg_env
export PYTHONPATH=/home/alon/Code/RMG-Py-i123c-audit:$PYTHONPATH
python -c "import rmgpy; print(rmgpy.__file__)"     # must be under this worktree
make build                                          # never bare `make`

python -m pytest -m "not functional and not database"          # 2 failed (pre-existing), 3020 passed

python docs/i123c-audit/probe_known_defects.py                 # exit 0
python docs/i123c-audit/probe_identity_consumers.py            # exit 0
python docs/i123-integration/probe_lithium_charge_network.py   # exit 0

mkdir -p /tmp/i123c-deck && cp docs/i123-integration/input.py /tmp/i123c-deck/
python rmg.py /tmp/i123c-deck/input.py \
  > >(tee /tmp/i123c-deck/stdout.log) 2> >(tee /tmp/i123c-deck/stderr.log >&2)   # RMG_EXIT=1, by design
sed -n '/^REACTIONS/,/^END/p' /tmp/i123c-deck/chemkin/chem.inp                   # both channels

python docs/i123c-audit/probe_readback_identity.py /tmp/i123c-deck               # exit 0
python docs/i123c-audit/probe_rewrite_units.py     /tmp/i123c-deck               # exit 1 -- THE BLOCKER
```

### 15.1 The scratch worktrees this pass created

Four, all detached, all disposable once this report is read. They exist because the pass had to
measure eight contributing branches and one baseline without writing to any tree it does not own:

| path | repository | at | why |
|---|---|---|---|
| `/home/alon/Code/RMG-Py-i123c-step` | RMG-Py | detached, last left at `948a783fd` | the stepping worktree for the per-branch sweep (§7) |
| `/home/alon/Code/RMG-Py-i123c-baseline` | RMG-Py | detached `a61dc1303` | the authoritative baseline, for the negative control (§12) and the read-refusal control (§8.4) |
| `/home/alon/Code/RMG-database-i123c-baseline` | RMG-database | detached `fb3c13c60` | the database baseline (§5.1) |
| `/home/alon/Code/RMG-database-i123c-step1` | RMG-database | detached `c7bd96292` | the first-pass database, to pair with `i126-chemkin-electrons` (§7) |

Each carries an untracked `rmgrc` written by this pass and nothing else; none is on a branch.
`git worktree remove --force <path>` on all four is safe.

