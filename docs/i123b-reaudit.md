# I-123, second pass — does the rebuilt union integrate?

**Verdict: NOT READY.** Two blocking failures, neither of them the one the first pass found.

1. **§7.1 — the radiative recombination sink never reaches the model.** RMG discards it as a
   duplicate of the electron-impact ionisation, because in the canonical representation the two
   are each other's reverse in heavy species and the duplicate check never looks at
   `Reaction.electrons`. The deck runs to completion and writes a mechanism with an ionisation
   **source and no sink**: the cation has no loss channel at all. This is the network the backlog
   exists to deliver, and the generated mechanism does not contain it.
2. **§7.2 — the emitted mechanism cannot be read back.** The Chemkin file RMG writes is refused by
   RMG's own Chemkin reader and by Cantera's `ck2yaml`, both on the `TDEP/` auxiliary line RMG's
   writer emits. Read-back is what the first pass could not test and what this pass was asked to
   prove. Removing that one line makes both readers accept the file, and the native Cantera YAML
   RMG writes alongside it reads back correctly — so the mechanism is not lost, but the Chemkin
   artifact does not round-trip.

The first pass's blocker — the Chemkin writer refusing every plasma rate law — **is fixed and
confirmed fixed** (§7.0). `rmgpy/chemkin.pyx` is byte-identical to `plasma`; the repair is entirely
in the shared converter, and the guard that used to refuse now passes.

Ten merges, five in each repository, all present and all conflict-free. Nothing was merged to
`plasma`, nothing was pushed, no PR was opened.

| | |
|---|---|
| RMG-Py | `/home/alon/Code/RMG-Py-i123b-reaudit`, branch `i123b-reaudit`, tip `9874fc706` |
| RMG-database | `/home/alon/Code/RMG-database-i123b-reaudit`, branch `i123b-reaudit-db`, tip `2d7123b08` |
| Shared branch (baseline) | RMG-Py `plasma` @ `a61dc1303` + RMG-database `plasma` @ `fb3c13c60` |
| Build | `make build` — **exit 0**; bare `make` **refused**, exit 2 |
| RMG-Py unit suite | baseline **2623 passed / 0 failed**, union **2834 passed / 2 failed** |
| RMG-database suite | baseline **41 passed / 4 failed**, union **173 passed / 4 failed** (same 4) |
| Union-only failures | **exactly one**, §5.2 |
| Contributing branches | 8 green, **2 red on their own branch**, §6 |
| Merge-clean? | **No** — two configuration commits must be reverted first, §2 |

Evidence labels: **[R]** code (file:line) · **[D]** database (file, entry) · **[M]** measured
(command and output shown) · **[I]** inference with its basis · **UNKNOWN** where nothing was
established. Raw logs referenced below live under `docs/i123b-reaudit/evidence/`.

---

## 1. The merge inventory, from git

### RMG-Py — `git log --oneline --graph plasma..HEAD` [M]

```
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

Verbatim, from `evidence/graph-rmgpy.log`. The range excludes `plasma` itself, whose tip
`a61dc1303 "Merge the reverse-reconstruction refusal (I-098)"` is the merge base of every line
above.

### RMG-database — `git log --oneline --graph plasma..HEAD` [M]

```
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

`evidence/graph-db.log`.

### Every merge's second parent is still its branch tip [M]

Checked rather than assumed — a branch that moved after being merged would make the audited union
differ from the branches the manager thinks were audited:

```
$ for m in f087cdf04 5e1686bcb 9f0eae40b c1f23e93e 582fd6ab9; do echo "$m p2=$(git rev-parse --short $m^2)"; done
f087cdf04  p2=99602b099   i102-quarantine        = 99602b099   MATCH
5e1686bcb  p2=ab3d84072   i112-marcus-work-terms = ab3d84072   MATCH
9f0eae40b  p2=ceb9fd5f1   i115-preflight-deck    = ceb9fd5f1   MATCH
c1f23e93e  p2=92e2d4234   i119-rr-registry       = 92e2d4234   MATCH
582fd6ab9  p2=fa042c08d   i110-make-guard        = fa042c08d   MATCH
```

Same in the database repository — `c7bd96292^2 = aabc3c622` (i119-recombination),
`1fb224371^2 = fcd87ce6e` (i114-ionisation-declaration), `ac944618e^2 = fc7bc101b`
(i104-alkali-plasma), `79320d363^2 = 069ea3d0b` (i103-electrochem-provenance),
`3b61056c1^2 = 1fcbb80fa` (i111-sei-reclassification). All ten match [M].

### One contribution that is not one of the ten [M]

`2d7123b08 "Give argon a sourced cation thermochemistry, and audit the gap it came from"` sits on
the database branch as a **direct commit, not a merge**. It is the tip of the separate branch
`i127-argon-thermo` (`git rev-parse i127-argon-thermo` = `2d7123b08`), so the pair under audit is
ten merges **plus** the export fix **plus** this. It adds `input/thermo/libraries/
PlasmaCationThermo.py` (170 new lines, one entry, `[Arp]`, `H298 = 1520.581 kJ/mol`) and four
database test files. It is included in every number below.

### Live work adjacent to this audit [M]

`git worktree list` shows `/home/alon/Code/RMG-Py-i126-chemkin-electrons` and
`/home/alon/Code/RMG-Py-i130-rms-plasma` both sitting at `1b79bee5b`, and
`/home/alon/Code/RMG-database-i127-argon-thermo` and `/home/alon/Code/RMG-database-i129-lithium-cation`
both at `2d7123b08` — that is, four other worktrees are based directly on this audit's two base
commits. Nothing in either was touched.

---

## 2. The pinned-configuration hazard: **still present, and now doubled**

**The branch is NOT merge-clean.** Two commits must come out before any merge to `plasma`.

| commit | what it does |
|---|---|
| `ae52de38a` | first pass: repointed `rmgrc` to `../RMG-database-i123-integration/input` |
| `9874fc706` | **this pass**: repointed it again, to `../RMG-database-i123b-reaudit/input` |

The first pass's commit is still present [M] — it was never reverted, and its target worktree is
the *previous* pass's database, one commit behind the pair under audit. Left as it was, a suite run
from this repo root would have loaded the wrong database, which is the same class of meaningless
green the first pass warned about. I repointed it rather than reverting it, because the audit needs
*a* pinned database and it must be this pass's; the revert requirement is unchanged and now covers
both commits.

The full delta is one line plus its comment [M]:

```
$ git diff plasma..HEAD -- rmgrc
-database.directory = ../RMG-database-plasma/input
+database.directory = ../RMG-database-i123b-reaudit/input
```

**The revert is a single command**, `git checkout plasma -- rmgrc`, and it is the only tracked
configuration difference between `plasma` and this branch — I checked the whole diff
(`git diff --stat plasma..HEAD`, 96 files) and no other tracked file is a harness pin [M].

If it is not reverted, `plasma` ships pointing at a worktree that exists only on this machine.

---

## 3. Which database actually loaded — measured, not trusted

Every probe in this report prints its resolved path at the head. The claim is never taken from
`rmgrc`:

```
$ python -c "from rmgpy import settings; import os; p=settings['database.directory']; \
             print(p); print(os.path.realpath(p)); print(os.path.isdir(p))"
/home/alon/Code/RMG-database-i123b-reaudit/input
/home/alon/Code/RMG-database-i123b-reaudit/input
True

$ python -c "import rmgpy, rmgpy.molecule.molecule as m; print(rmgpy.__file__); print(m.__file__)"
/home/alon/Code/RMG-Py-i123b-reaudit/rmgpy/__init__.py
/home/alon/Code/RMG-Py-i123b-reaudit/rmgpy/molecule/molecule.cpython-39-x86_64-linux-gnu.so
```

[M] — and the compiled `.so`, not a stray `.py`, so the measurements are of the built tree.
The same two lines appear at the head of `evidence/probe-lithium-charge-network.log`,
`evidence/probe-sink-collapse.log`, `evidence/probe-roundtrip.log`, `evidence/rms-loud-failure.log`,
and in every per-branch summary under `evidence/branch-*.summary.log`.

The baseline was resolved the same way: `/home/alon/Code/RMG-Py-i123b-baseline/rmgpy/__init__.py`
with `database.directory = /home/alon/Code/RMG-database-plasma/input` [M], which is the shared
branch paired with itself. `/home/alon/Code/RMG-database-plasma` was confirmed at `fb3c13c60`, the
`plasma` tip, before being used, and only read from.

---

## 4. Build, before any test result

### `make build` — exit 0 [M]

```
$ make build
... (72 KB of cython + gcc, evidence/make-build-union.tail.log)
MAKE_BUILD_EXIT=0
```

104 extension modules built in the union worktree.

### Bare `make` — still refused, exit 2 [M]

```
$ make
Refusing to modify the shared RMG environment.
Use `make build` for an in-place worktree build.
Editable installation requires an explicit maintenance procedure:
    make unsafe-install-shared-env CONFIRM_SHARED_ENV_MUTATION=yes
make: *** [Makefile:72: guard] Error 1
$ echo $?
2
```

`make install` is refused identically, exit 2 [M] (`evidence/bare-make.log`,
`evidence/bare-make-install.log`). The guard works.

### An unadvertised build difference between the two branches [M]

On `plasma`, `make build` is `python setup.py build_ext --inplace` with **no dependency on
`check-pydas`** (`Makefile:38`), so in a fresh worktree it dies:

```
Cython.Compiler.Errors.InternalError: Internal compiler error: 'settings.pxi' not found
make: *** [Makefile:39: build] Error 1
```

The union's `i110-make-guard` merge is what adds that dependency. Every baseline and
contributing-branch measurement below therefore ran `python utilities.py check-pydas` explicitly
first; that is a deviation from "same command", and it applies only to the build step, not to the
test step. This is a small argument *in favour* of the merge, not against it.

---

## 5. The unit suites

Same command on both, verbatim, the one `make test` runs:

```
python -m pytest -m "not functional and not database"
```

| tree | passed | failed | skipped | deselected | wall |
|---|---|---|---|---|---|
| shared branch (`plasma` + database `plasma`) | **2623** | **0** | 50 | 83 | 8:53 |
| union (`i123b-reaudit` + `i123b-reaudit-db`) | **2834** | **2** | 50 | 131 | 9:53 |

[M] — `evidence/unit-suite-baseline.summary.log`, `evidence/unit-suite-union.summary.log`. The union
adds 211 passing tests and two failures.

Database suites, `python -m pytest test -q`:

| tree | passed | failed |
|---|---|---|
| database `plasma` + RMG-Py `plasma` | **41** | **4** |
| database union + RMG-Py union | **173** | **4** |

[M] — `evidence/db-suite-baseline-pinned.log`, `evidence/db-suite-union.log`. The four are the same
four in both, all in `test_plasma_electron_attachment.py`:
`test_trained_species_resolve_to_their_own_rate_through_their_own_node[O2|OH|O]` and
`test_o2_rate_comes_from_the_training_set_not_a_library`. **Pre-existing on the shared branch, not
caused by the merge, and not diagnosed here** — they were outside this pass's Verifier.

> A finding on the way to that baseline: on the `plasma` database branch the suite has **no
> `database.directory` pin**, so it resolves to `../RMG-database/input` — the polymer database — and
> **45 of 45 tests error out** before running (`evidence/db-suite-baseline-unpinned.log`). The
> baseline row above was obtained with an untracked scratch `rmgrc` pinning the worktree. The union
> fixes this properly, in `bdbbf2864 "test: pin database.directory for the suites in this
> repository"` (`test/conftest.py`, 29 lines). Another argument for the merge.

### 5.1 The two union failures, named

```
FAILED test/rmgpy/electronPlacementTest.py::TestElectronPlacementResolver::
       test_declaration_registry_is_explicit_and_closed
  assert {'Plasma_Electron_Attachment': (1,0), 'Cation_R_Recombination': (1,0),
          'PlasmaRadiativeRecombination': (1,0), 'PlasmaElectronImpactIonization': (1,2)}
      == {'Plasma_Electron_Attachment': (1,0), 'Cation_R_Recombination': (1,0),
          'PlasmaElectronImpactIonization': (1,2)}

FAILED test/rmgpy/preflightDeckFamilyExclusionTest.py::PlasmaDeckFamilyExclusionTest::
       test_plasma_deck_does_not_declare_family[docs/i102-quarantine/input.py]
  plasma deck docs/i102-quarantine/input.py declares Cation_R_Recombination, a lithium-ion-battery
  SEI family excluded from every plasma configuration.
```

### 5.2 Union-only failures: **exactly one**, and it is still there

The first pass found exactly one genuine union-only failure. So does this pass, and it is the same
one.

**Union-only — `preflightDeckFamilyExclusionTest[docs/i102-quarantine/input.py]`.** Attributed to
the pair **i115-preflight-deck × i102-quarantine** [M]:

```
$ for b in ...; do git cat-file -e $b:test/rmgpy/preflightDeckFamilyExclusionTest.py; \
                  git cat-file -e $b:docs/i102-quarantine/input.py; done
  i115-preflight-deck   : test=present  deck=absent
  i102-quarantine       : test=absent   deck=present
  i119-rr-registry      : test=absent   deck=absent
  i112-marcus-work-terms: test=absent   deck=absent
  i110-make-guard       : test=absent   deck=absent
  plasma                : test=absent   deck=absent
```

The test discovers decks by globbing the tree (`preflightDeckFamilyExclusionTest.py:165`,
`ALL_DECKS = _discover_decks()`) rather than from a fixed list [R], so it can only see i102's deck
once both branches are in one tree. Neither branch is wrong on its own; the sweep from one meets a
deck from the other. Mechanically: i102's reproduction deck declares `Cation_R_Recombination`, which
i115 established must not appear in a plasma configuration.

**NOT union-only — `test_declaration_registry_is_explicit_and_closed`.** This one is **inherited**:
it already fails on `i119-rr-registry` by itself [M]:

```
$ cd /home/alon/Code/RMG-Py-i123b-b-i119 && python -m pytest test/rmgpy/electronPlacementTest.py -q
FAILED test/rmgpy/electronPlacementTest.py::TestElectronPlacementResolver::
       test_declaration_registry_is_explicit_and_closed
1 failed, 50 passed, 2 warnings in 1.38s
```

Within that one seven-ticket branch, `92e2d4234 "Let a kinetics library own radiative recombination"`
adds a fourth registry entry and the closure test's expected set was left at three. It is the
branch's own defect, carried into the union unchanged. It is a stale assertion, not a behaviour
change — the resolver's own behaviour tests all pass, and the live registry with four entries is what
the working probe in §7 exercises.

---

## 6. Each contributing branch's own green/red state, measured this pass

Every RMG-Py branch was checked out detached into its own worktree, `check-pydas`'d, `make build`'d,
and run with the same suite command. **Eight of the ten are green; two are red.**

### RMG-Py [M] — `evidence/branch-*.summary.log`

| branch | tip | resolved database | result |
|---|---|---|---|
| `i110-make-guard` | `fa042c08d` | `/home/alon/Code/RMG-database/input` | **RED** — 2 failed, 23 errors |
| `i119-rr-registry` | `92e2d4234` | `RMG-database-plasma/input` | **RED** — 1 failed, 2684 passed |
| `i115-preflight-deck` | `ceb9fd5f1` | `RMG-database-plasma/input` | GREEN — 2630 passed |
| `i112-marcus-work-terms` | `ab3d84072` | `RMG-database-plasma/input` | GREEN — 2667 passed |
| `i102-quarantine` | `99602b099` | `RMG-database-plasma/input` | GREEN — 2641 passed |

**`i110-make-guard` is red for a configuration reason, not a code reason.** It carries **no tracked
`rmgrc`** — it is the branch that changes the Makefile, and it predates the pin — so RMG falls back
to `../RMG-database/input`, the **polymer** database. All 23 errors are one message,
`DatabaseError: Unrecognized item "Plasma_Electron_Attachment"`, and the two failures are
`thermoTest.py::test_parse_thermo_comments` (`KeyError: 'C=XR2'`) and
`test_adsorbate_thermo_raises_error` (`KeyError: 'NX'`) — polymer-database group names. Pinned to the
plasma database it is **fully green** [M]:

```
$ printf 'database.directory = /home/alon/Code/RMG-database-plasma/input\n' > rmgrc
$ python -m pytest -m "not functional and not database"
=== 2635 passed, 50 skipped, 68 deselected, 54 warnings in 485.67s (0:08:05) ===
```

This is exactly the hazard of §2 seen from the other side, and it is harmless after the merge: the
union has an `rmgrc`.

**`i119-rr-registry` is genuinely red**, on the stale closure assertion of §5.2. This is the branch
the first pass found red; it is still red, and nothing has been done about it.

> **A measurement error of my own, corrected.** I ran the five branch suites in parallel. Three
> reported an extra failure, `qmMainTest.py::test_get_thermo_data_mopac`, with
> `FileNotFoundError: /home/alon/Code/testing/qm/QMscratch/...` — a **fixed scratch path outside
> every worktree**, which five concurrent runs share and clobber. Re-run serially, that test passes
> on all five [M] (`evidence/qm-serial-recheck.log`), and the full suites of the two affected
> branches were re-run serially before the table above was written
> (`evidence/branch-i102-serial.summary.log`, `evidence/branch-i110-pinned.summary.log`). The table
> reports the serial numbers. Two consequences worth recording: **RMG's QM test suite is not safe to
> run concurrently across worktrees**, and any past "green" or "red" obtained from a parallel sweep
> is suspect.

### RMG-database [M] — `evidence/db-contributing-branches.log`

A database branch has no green/red state on its own: it needs an RMG-Py to be measured against, and
none of these branches records which. Measured against both candidates:

| branch | tip | vs RMG-Py union | vs RMG-Py `plasma` |
|---|---|---|---|
| `i111-sei-reclassification` | `1fcbb80fa` | 4 failed, 59 passed | 4 failed, 59 passed |
| `i103-electrochem-provenance` | `069ea3d0b` | 4 failed, 41 passed | 4 failed, 41 passed |
| `i104-alkali-plasma` | `fc7bc101b` | 4 failed, 41 passed | 4 failed, 41 passed |
| `i114-ionisation-declaration` | `fcd87ce6e` | 4 failed, 71 passed | 13 failed, 62 passed |
| `i119-recombination` | `aabc3c622` | 4 failed, 107 passed | 24 failed, 87 passed |

All five carry exactly the four pre-existing attachment failures and nothing else, **provided they
are paired with the RMG-Py union**. Against the shared branch's RMG-Py, `i114` and `i119` add 9 and
20 further failures — they are the database halves of RMG-Py feature work and do not stand alone.
That is expected, not a defect, but it does mean **neither can be merged to the database `plasma`
before its RMG-Py counterpart lands**.

So: green, on the pairing that matters, all five.

---

## 7. The lithium charge network, end to end

### 7.0 First — the first pass's blocker is gone [M]

`rmgpy/chemkin.pyx` is **byte-identical** between `plasma` and the union
(`git diff plasma..HEAD -- rmgpy/chemkin.pyx` is empty) [M]. The writer and its refusal
(`chemkin.pyx:2060`, `check_electron_reactant_order`) are untouched. The fix is entirely in
`rmgpy/electron_balance.py` and `rmgpy/electron_placement.py`, and it works: the deck's `chem.inp`
now carries

```
Li(2)+e-(1)=>[Lip](3)+e-(1)+e-(1)                   1.490e+18 -0.267    162.150
    TDEP/e-(1)/   ! VoronovEIArrhenius exported as a modified-Arrhenius fit of k(Te) ...
```

— one electron in, two out, second order, which is what the Voronov coefficient's
cm³ molecule⁻¹ s⁻¹ dimensionality requires. The guard that used to kill the run passes it. **The
first pass's §6.1 is closed.**

### 7.1 Stages 1–4: load, balance, placement, reactor — all pass [M]

`python docs/i123-integration/probe_lithium_charge_network.py`, full output in
`evidence/probe-lithium-charge-network.log`, exit 0:

| stage | result |
|---|---|
| **1 LOAD** | both libraries load out of the union database; one entry each; the source produces exactly what the sink consumes |
| **2 BALANCE** | `is_balanced` true for both; each closes its own charge gap |
| **3 PLACEMENT** | both are declared owners; ionisation resolves to `Li + e⁻ ⇒ Li⁺ + 2 e⁻` (1,2); recombination to `Li⁺ + e⁻ ⇒ Li + hν` (1,0); both second order; the canonical reactions are not mutated |
| **4 REACTOR** | one `PlasmaReactor` accepts both; `kf[source] = 5.031584e+07` matches Voronov, `kf[sink] = 1.468188e+05` matches Badnell; the cation has a non-zero loss channel |

The machinery is sound. **Everything after this point is about what RMG does with it.**

### 7.2 Stage 5, the deck: it runs, and the sink is not in the answer — **BLOCKING**

```
$ python rmg.py docs/i123-integration/input.py
RMG_EXIT=0
The final model core has 7 species and 1 reactions
The final model edge has 0 species and 0 reactions
```

One reaction. The `chem.inp` contains the ionisation and nothing else; `chem_edge.inp` is
byte-identical to it. `RMG.log` records the reason in a single informational line [M]:

```
Adding reaction library PlasmaRadiativeRecombination to model edge...
This library reaction was not new: [Lip](3) => Li(2)
```

The sink was **discarded as a duplicate of the source**. I probed this against the running code
rather than reading it — `python docs/i123b-reaudit/probe_sink_collapse.py`,
`evidence/probe-sink-collapse.log`, exit 0 [M]:

| stage | result |
|---|---|
| 1 | the two are physically distinct: `electrons = +1` vs `−1`; their heavy-species halves are exact mirrors, `[Li] ⇒ [Lip]` vs `[Lip] ⇒ [Li]` |
| 2 | source added first ⇒ **sink** judged not new, `is_new=False`, dropped |
| 3 | sink added first ⇒ **source** judged not new, dropped. **Order-dependent**: whichever is added second disappears |
| 4 | forcing `Reaction.electrons` to 0, −1 or 12345 changes nothing — the verdict never consults it |

Stage 3 is what rules out the benign reading. This is not RMG judging the sink unimportant; it is
RMG unable to tell the two reactions apart, and keeping whichever it saw first.

The mechanism is at `rmgpy/rmg/model.py:504-522` [R] — the cross-library duplicate check matches
`rxn_id == rxn_id0[::-1]`, i.e. the reverse direction, and `grep electron` over
`check_for_existing_reaction` and `make_new_reaction` (lines 424–600) returns nothing [M]. In the
canonical representation the electron is the scalar `Reaction.electrons` and is not a participant,
so ionisation and radiative recombination are literally each other's reverse in every field the
check looks at.

**Why this blocks.** The deliverable is a lithium plasma mechanism with a cation source *and* a
cation loss channel. What RMG writes has a source and no loss channel: Li⁺ accumulates
monotonically, and any downstream simulation of the written mechanism is wrong in the direction that
matters. The probe in §7.1 proves the reactor *can* hold both — but the reactor is fed by the model,
and the model has already dropped one. The failure is silent at the exit code (`RMG_EXIT=0`) and
carries no warning; the only trace is an `INFO` line that reads like routine deduplication.

This defect is reachable only on the union — the two libraries are never both present anywhere else
— but it is not an artifact of merging: it is a property of RMG's duplicate check meeting the
canonical electron representation, and it would appear on any tree carrying both libraries.

### 7.3 Stage 6, read-back — **BLOCKING**

`python docs/i123b-reaudit/probe_roundtrip.py <deck>`, `evidence/probe-roundtrip.log` [M]:

| trip | reader | outcome |
|---|---|---|
| A | RMG's `load_chemkin_file` on the file RMG just wrote | **REFUSED** — `ChemkinError: 'e-(1)' doesn't look like a collision efficiency for species TDEP in line 'TDEP/E-(1)/'` (`chemkin.pyx:632`) |
| B | same file, `TDEP/` line removed | **READ** — 7 species, 1 reaction, `Li(2) + e-(1) => [Lip](3) + e-(1) + e-(1)` |
| C | Cantera on RMG's own `cantera2/chem_annotated.yaml` | **READ** — 7 species, 1 reaction, same equation |

And Cantera's Chemkin importer agrees with RMG's [M] (`evidence/ck2yaml-tdep.log`):

```
--- ck2yaml on chem.inp
  REFUSED: InputError -- could not convert string to float: 'e-(1)'
--- ck2yaml on chem_no_tdep.inp
  CONVERTED: 7 species, 1 reactions
     Li(2) + e-(1) => [Lip](3) + e-(1) + e-(1)
```

RMG hits this **inside its own run**: the final "Translating final chemkin file into Cantera yaml"
step fails, `cantera_from_ck/` is left empty, and the traceback is written to `RMG.log` as
`Error: Could not generate Cantera files for some reason` — while `rmg.py` still exits **0** [M].

**Why this blocks.** Verifier item 8 required write-out *and read-back*, and read-back of the
Chemkin artifact does not work. Two honest qualifications, both decision-relevant: the reader is
**unchanged by this merge** (`chemkin.pyx` byte-identical to `plasma`), so the merge does not create
the reader gap — it creates the first files that hit it; and the mechanism itself is **not lost**,
because the native Cantera YAML round-trips correctly and the Chemkin file round-trips as soon as
one line is removed. The electron semantics do degrade on the way back — the reloaded reaction is a
plain `Arrhenius` with `electrons=0` — which is the documented lossiness of the Chemkin form, not a
new defect.

---

## 8. The two known non-blocking defects — confirmed as described

### 8.1 The RMS YAML writer refuses these rate laws, **loudly** [M][R]

The property the decision rests on is that it fails at save rather than emitting wrong numbers.
Confirmed at both levels.

At the writer (`evidence/rms-loud-failure.log`):

```
--- PlasmaElectronImpactIonization -> [Li] => [Lip]   kinetics=VoronovEIArrhenius
  LOUD: raised ValueError
  message: Object of type <class 'rmgpy.kinetics.arrhenius.VoronovEIArrhenius'> does not have a
           defined conversion to ReactionMechanismSimulator format
  at: File ".../rmgpy/yaml_rms.py", line 252, in obj_to_dict
--- PlasmaRadiativeRecombination -> [Lip] => [Li]   kinetics=BadnellRRArrhenius
  LOUD: raised ValueError  (same message, BadnellRRArrhenius)
```

`yaml_rms.py:251-253` is a terminal `else: raise ValueError(...)` in the type dispatch [R] — there is
no fall-through that could emit a default, so a silent wrong number is structurally impossible here.

At the deck, with `generateRMSYAML=True` restored (`evidence/deck-rms-enabled.stderr.tail.log`):

```
$ python rmg.py <deck with generateRMSYAML=True>
RMG_RMS_EXIT=1
  File ".../rmgpy/rmg/main.py", line 2254, in save_everything
  File ".../rmgpy/yaml_rms.py", line 252, in obj_to_dict
ValueError: Object of type <class '...VoronovEIArrhenius'> does not have a defined conversion to
            ReactionMechanismSimulator format
```

Uncaught, exit 1, dies in `save_everything()` — the same site the first pass's Chemkin blocker died
in. **Loud, not silent.** The deck-level workaround is in place and documented in the deck itself
(`docs/i123-integration/input.py:123-132`, `generateRMSYAML=False` with the reason written beside
it), and with it the run exits 0. Acceptable as briefed. Not touched.

### 8.2 The lithium cation enthalpy is wrong, and **identically wrong on the shared branch** [D][M][I]

`input/thermo/libraries/LithiumPrimaryThermo.py`, entry `index = 65`, `label = "[Lip]"`:
`E0 = (526.738, 'kJ/mol')` [D]. This is the entry the deck actually uses — the emitted
`chem_annotated.inp` attributes `[Lip](3)` to `! Thermo library: LithiumPrimaryThermo` [M].

The gap: ΔfH°(Li⁺, g) at 0 K is ΔfH°(Li, g) + IE(Li) ≈ 159.3 + 520.2 ≈ 679.5 kJ/mol, against the
entry's 526.7 — short by ≈ 153 kJ/mol ≈ **1.58 eV**, consistent with the ~1.5 eV briefed [I, basis:
standard atomisation enthalpy and first ionisation energy of lithium; not re-derived from a source
here]. The shortfall is close to the atomisation enthalpy itself, consistent with the entry being
referenced to Li(g) rather than Li(s) — noted, not diagnosed.

**The merge neither creates nor worsens it** [M]:

```
$ git diff --stat plasma..HEAD -- input/thermo/libraries/LithiumPrimaryThermo.py
(no output — the file is byte-identical)
```

Not corrected here, as instructed. Its ticket owns it.

---

## 9. Negative control — clean [M]

`examples/rmg/minimal/input.py`, byte-identical on both trees, run on each with its own pair:

```
union    : RMG-Py-i123b-reaudit  + RMG-database-i123b-reaudit  -> exit 0
baseline : RMG-Py-i123b-baseline + RMG-database-plasma         -> exit 0
The final model core has 26 species and 66 reactions      (both)
The final model edge has 146 species and 332 reactions    (both)

  IDENTICAL  chem.inp                 md5 308924af2ca69954eddcba1853080886
  IDENTICAL  chem_edge.inp            md5 6fc758854b2a3975c0c3d5216f7e8523
  IDENTICAL  chem_annotated.inp       md5 d5108f38ff4634aec225776ae82c47ff
  IDENTICAL  species_dictionary.txt   md5 c321b3b67de7ad4ac3bddfb23beea391
  IDENTICAL  tran.dat                 md5 d78237511ff04ffef29d039f66d97eb8
```

`evidence/negative-control.log`. Including the annotated file, which carries thermo and kinetics
provenance comments and would show any change of source. An ordinary non-plasma mechanism builds,
runs and exports byte-identically across the merge.

Consistent with the database delta, which touches only `Cation_R_Recombination`, the two new plasma
kinetics libraries, the new `PlasmaCationThermo.py`, docs and tests [M].

---

## 10. Verdict

**NOT READY.**

**Blocking failure 1 — the generated mechanism has no cation sink (§7.2).** RMG discards the
radiative recombination as a duplicate of the ionisation, silently, order-dependently, and without
ever consulting the field that distinguishes them. The lithium plasma deck runs to completion and
writes a one-reaction mechanism in which Li⁺ has no loss channel. This is the deliverable, and it is
not delivered. `rmgpy/rmg/model.py:504-522`.

**Blocking failure 2 — the emitted mechanism does not round-trip (§7.3).** Neither RMG's Chemkin
reader nor Cantera's `ck2yaml` can read the file RMG writes, both on the `TDEP/` line; RMG's own
final conversion step fails inside the run and leaves `cantera_from_ck/` empty while still exiting
0. Read-back was an explicit requirement of this pass. Mitigating, and worth weighing: the reader is
byte-identical to `plasma`, the native Cantera YAML reads back correctly, and removing one line
makes both Chemkin readers accept the file.

**Not blocking, and not re-litigated:** the RMS writer gap (§8.1, confirmed loud) and the Li⁺
enthalpy (§8.2, confirmed identical to the shared branch).

**Also required before merge, independent of the above:** revert both configuration commits, §2.

**What the merge does deliver**, so the cost of *not* merging is on the record: +211 passing unit
tests and +132 passing database tests against the shared branch; a working `make` guard; a database
suite that can find its own database; the export converter fix that closes the first pass's blocker;
and a lithium charge network that loads, balances, resolves and evaluates correctly in the reactor
(§7.1). The two blockers are both downstream of that, in the model builder and in the file readers.

---

## 11. What this pass could not reach

- **Whether the §7.2 sink collapse also affects family-generated plasma reactions**, as opposed to
  the library-to-library case measured here. Not probed — but `rmgpy/rmg/model.py:485` applies the
  same `rxn_id == rxn_id0[::-1]` reverse match to `KineticsFamily` reactions [R], so the shape is
  there. Worth one probe before anyone assumes the blast radius is one library pair.
- **A fix, or a design for one, for either blocker.** Out of scope by the brief; the report names
  the mechanism and stops.
- **The four pre-existing `Plasma_Electron_Attachment` database failures.** Present identically on
  the shared branch; not diagnosed. They were the first pass's open item and remain open.
- **The `functional` and `database`-marked pytest suites** (131 and 83 deselected respectively) and
  **`test/regression/`**. Only `-m "not functional and not database"` was run, on every tree, which
  is what "the full unit suite" means here — but it is not everything that CI runs.
- **Any platform but this one.** Linux, Python 3.9, one machine. No macOS, no other Python.
- **The RMS runtime path.** Julia was never invoked; §8.1 measures the writer only.
- **Whether the deck's transport and diffusion numbers are right.** Untouched; that is I-121/I-128
  territory and nothing here bears on it.
- **The reason `i119-rr-registry` was left red** after the first pass reported it. Measured again,
  not investigated.
- **Whether `plasma` is the right merge target ref at all.** Local `plasma` is 10 commits *ahead* of
  `origin/plasma` [M]; the first pass baselined against `origin/plasma` and this one against local
  `plasma`. Which is authoritative was the first pass's open question and is still unanswered.

---

## 12. Reproducing this

```bash
conda activate rmg_env
cd /home/alon/Code/RMG-Py-i123b-reaudit
export PYTHONPATH=/home/alon/Code/RMG-Py-i123b-reaudit:$PYTHONPATH
python utilities.py check-pydas && make build

python docs/i123-integration/probe_lithium_charge_network.py   # §7.1, exit 0
python docs/i123b-reaudit/probe_sink_collapse.py               # §7.2, exit 0
python rmg.py docs/i123-integration/input.py                   # the deck
python docs/i123b-reaudit/probe_roundtrip.py <deck-output-dir> # §7.3
python -m pytest -m "not functional and not database"          # §5
```

Scratch worktrees used for the baseline and the ten contributing branches
(`RMG-Py-i123b-baseline`, `RMG-Py-i123b-step`, `RMG-Py-i123b-b-*`,
`RMG-database-i123b-baseline`, `RMG-database-i123b-step`) were created detached, never committed
to, and removed afterwards. No branch ref was moved. No shared checkout was written to.
