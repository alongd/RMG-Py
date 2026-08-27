# I-123, fourth pass — the union with all three previous blockers cleared

**Verdict: NOT READY.** One blocking failure, in a path no previous pass exercised: **the seed
mechanism RMG writes for every plasma run cannot be restarted from.** The generated
`restart_from_seed.py` builds a model with the ionisation channel counted **twice** and then dies
at the first Chemkin save, exit 1. An ordinary non-plasma mechanism restarts cleanly through the
same code, so this is specific to the electron representation this union delivers. §9.

Everything the brief named as the deliverable **does work**, measured here and not inherited: the
lithium charge network runs end to end, three successive Chemkin round trips move no rate constant
by a single digit, `scripts/mergeModels.py` keeps both channels with both rate constants exact, and
an ordinary mechanism exports byte-identically to the shared branch. The three previous blockers
are all genuinely repaired on this base. §§7, 8, 11.

Evidence labelling: **[R]** code (file:line) · **[D]** database (file, entry) · **[M]** measured
(command and output shown) · **[I]** inference with its basis.

Raw output for every measurement is under [`docs/i123d-audit/evidence/`](i123d-audit/evidence/);
the three probes this pass wrote are [`probe_seed_roundtrip.py`](i123d-audit/probe_seed_roundtrip.py),
[`probe_merge_models.py`](i123d-audit/probe_merge_models.py) and
[`branch_sweep.sh`](i123d-audit/branch_sweep.sh).

---

## 0. Two things wrong with the brief's own premises, found before any measurement

Reported rather than worked around, because both would have silently changed what this audit
measured.

**0.1 A required-reading file does not exist.** The brief names
`docs/contracts/closed/i145-units-laundering.md` as one of five documents to read on my own base.
There is no such file, and no `docs/contracts/` directory of any kind on this branch **[M]**:

```
$ ls docs/contracts/closed/i145-units-laundering.md
ls: cannot access 'docs/contracts/closed/i145-units-laundering.md': No such file or directory
$ find . -name '*i145*'
./docs/i145-units-laundering
```

The units work is documented as `docs/i145-units-laundering/` (two runnable probes plus evidence)
and in the three commit messages `8097027b1`, `56cb278fa`, `e5e42ae41`, which are unusually full.
I read those instead. Nothing was lost, but a reader of the brief would conclude a design document
exists that does not.

**0.2 The tracked `rmgrc` on my base pointed at another pass's database worktree.** Inherited on
`i123d-audit`, `database.directory` read `../RMG-database-i123c-audit/input` — the **third** pass's
worktree, which this pass does not own and which another branch is free to move under a run in
progress. Every probe and every suite here would have read its kinetics libraries out of a checkout
whose contents this pass cannot vouch for. Repointed at this pass's own sibling in `cd6cde425`.
This is the exact hazard §3 exists for, one pass later and one worktree along. See §3 for the
merge-clean consequence.

---

## 1. The authoritative base, and the local-versus-remote gap

The authoritative base is the **local** shared branch, not the remote. The gap is **unchanged**
from the third pass, in both repositories **[M]**:

```
$ git rev-parse --short plasma origin/plasma
a61dc1303
6d3c03823
$ git rev-list --count origin/plasma..plasma      # local ahead of remote
10
$ git rev-list --count plasma..origin/plasma      # remote ahead of local
0
```

Same two SHAs the third pass recorded, so the same ten commits: the I-098 reverse-reconstruction
refusal and the walls A–D chain. The third pass established that the correction moves no verdict
(its §1.4); nothing here re-opens that, and the number that would have carried the risk — the
baseline unit-suite count — comes out at **2623 passed / 0 failed**, identical to what passes 2 and
3 measured against the same local ref (§5).

The RMG-database repository has no divergence at all **[M]**:

```
$ git -C /home/alon/Code/RMG-database-i123d-audit rev-list --count origin/plasma..plasma ; \
  git -C /home/alon/Code/RMG-database-i123d-audit rev-list --count plasma..origin/plasma
0
0
```

Both are `fb3c13c60`. **Verdict-neutral, still.**

---

## 2. The merge inventory of each branch, from git

### 2.1 RMG-Py, branch `i123d-audit` @ `cd6cde425`

Full graph: [`evidence/graph-rmgpy.log`](i123d-audit/evidence/graph-rmgpy.log). **[M]**
`git log --oneline --graph plasma..HEAD`:

```
* cd6cde425 Point the fourth-pass audit branch's rmgrc at its own database worktree
* e5e42ae41 Record the units-pairing evidence: probes, before, after, red check
* 56cb278fa Stop VoronovEIArrhenius relabelling a metre-based A as centimetre-based
* 8097027b1 Pair the plasma Chemkin reduction's A-factor with its own units
* 7f1392214 Report the I-123 third audit: NOT READY, one blocking failure
* 34de41c76 Point the third-pass audit branch's rmgrc at its own database worktree
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
*   5e1686bcb Merge i112-marcus-work-terms: let the Marcus work terms reach the barrier they belong in
*   9f0eae40b Merge i115-preflight-deck: exclude the battery-SEI family from the plasma preflight decks
*   c1f23e93e Merge i119-rr-registry: the charge-network chain, seven tickets in one branch
*   582fd6ab9 Merge i110-make-guard: refuse a bare `make` that installs into the shared conda environment
```

Nine contributing branches, and **all nine tips are ancestors of HEAD** — nothing has drifted out
from under the merges since they were made **[M]**:

| branch | tip | in HEAD? |
|---|---|---|
| `i110-make-guard` | `fa042c08d` | ancestor |
| `i119-rr-registry` | `92e2d4234` | ancestor |
| `i115-preflight-deck` | `ceb9fd5f1` | ancestor |
| `i112-marcus-work-terms` | `ab3d84072` | ancestor |
| `i102-quarantine` | `99602b099` | ancestor |
| `i126-chemkin-electrons` | `1b79bee5b` | ancestor |
| `i134-duplicate-drops-sink` | `9b8cc38e7` | ancestor |
| `i135-tdep-roundtrip` | `948a783fd` | ancestor |
| `i145-units-laundering` | `e5e42ae41` | ancestor |

Note the shape: `i126`, and the whole of `i145`, are **direct commits on the union's first-parent
chain**, not merges. `i123d-audit` was created at the tip of `i145-units-laundering`, so before my
`rmgrc` commit the two branch names pointed at the same object. The fourth pass therefore adds no
new merge — it audits the third pass's union plus the units repair.

### 2.2 RMG-database, branch `i123d-audit-db` @ `d7c8595ce`

Full graph: [`evidence/graph-rmgdatabase.log`](i123d-audit/evidence/graph-rmgdatabase.log).

**The database half is byte-for-byte the third pass's** **[M]**:

```
$ git log --oneline i123c-audit-db..HEAD ; git log --oneline HEAD..i123c-audit-db
(both empty)
```

Twenty commits ahead of `plasma` @ `fb3c13c60`, over five contributing branches, all ancestors
**[M]** — `i103-electrochem-provenance` `069ea3d0b`, `i104-alkali-plasma` `fc7bc101b`,
`i111-sei-reclassification` `1fcbb80fa`, `i114-ionisation-declaration` `fcd87ce6e`,
`i119-recombination` `aabc3c622`, plus the direct commit `2d7123b08` (argon cation thermo) and the
third pass's report.

**Two database branches exist and are deliberately NOT in this union** **[M]**:
`i129-lithium-cation-enthalpy` @ `75e0cb6f0` and `i120-argon-recombination` @ `5f52259f6`, neither
an ancestor. The first is why §10.4's enthalpy gap is still open on this tree.

---

## 3. The pinned-configuration hazard — the branch is **NOT merge-clean**

`rmgrc` is tracked, and `database.directory` in it decides which database every run in this
worktree reads. On `plasma` the line must read `../RMG-database-plasma/input`. **Four** commits now
move it, one per audit pass **[M]**:

```
$ git log --oneline -- rmgrc
cd6cde425 Point the fourth-pass audit branch's rmgrc at its own database worktree      <- this pass
34de41c76 Point the third-pass audit branch's rmgrc at its own database worktree
9874fc706 Point the re-audit branch's rmgrc at the re-audit database worktree
ae52de38a Point the integration branch's rmgrc at the I-123 integration database
```

**Plainly: this branch is not merge-clean. All four of those commits must be reverted — or the one
line restored to `../RMG-database-plasma/input` — before `i123d-audit` goes anywhere near
`plasma`.** Merging as-is would repoint every plasma run on the shared branch at a private audit
worktree that will be deleted.

Why the pin has to be there while the audit runs: without it `database.directory` resolves to the
sibling `RMG-database`, which is the polymer-branch database and carries neither of the two
kinetics libraries this union exists to deliver. Every plasma suite would then fail with
`DatabaseError: Unrecognized item "Plasma_Electron_Attachment"`, which looks exactly like a
regression and is not one; worse, a green run against it would be meaningless.

**Configuration is not trusted anywhere in this report.** Every probe prints the resolved path at
the head of its own run, and every one of them printed **[M]**:

```
database.directory = /home/alon/Code/RMG-database-i123d-audit/input
resolved realpath  = /home/alon/Code/RMG-database-i123d-audit/input
rmgpy              = /home/alon/Code/RMG-Py-i123d-audit/rmgpy/__init__.py
compiled           = /home/alon/Code/RMG-Py-i123d-audit/rmgpy/molecule/molecule.cpython-39-x86_64-linux-gnu.so
```

The deck run's own header agrees, independently of the probes
([`evidence/deck-stdout.tail.log`](i123d-audit/evidence/deck-stdout.tail.log)):
`database.directory = /home/alon/Code/RMG-database-i123d-audit/input (from rmgrc)`, RMG-Py HEAD
`cd6cde425`, RMG-database HEAD `d7c8595ce`.

---

## 4. Build, stated before any test result

**[M]** `make build`, from a tree with **zero** `.so` files (`find . -name '*.so' | wc -l` → `0`):

```
### make build exit: 0 ###
real    3m38.544s
```

[`evidence/make-build.tail.log`](i123d-audit/evidence/make-build.tail.log). No errors on stderr
(`grep -ciE '\berror\b'` → `0`).

**The bare-`make` guard is still armed** **[M]**
([`evidence/bare-make-guard.log`](i123d-audit/evidence/bare-make-guard.log)):

```
$ make
Refusing to modify the shared RMG environment.
Use `make build` for an in-place worktree build.
...
make: *** [Makefile:72: guard] Error 1
### bare make exit: 2 ###
```

The baseline `plasma` tree also built cleanly (`### baseline make build exit: 0 ###`), and every
one of the ten trees in the §7 sweep reported `make build exit: 0`.

---

## 5. The unit suites — same command, run serially

`python -m pytest -m "not functional and not database"`, one at a time. Serial on purpose:
concurrent pytest runs across worktrees share a fixed scratch path outside every worktree and
produce false failures, which is how the second pass acquired phantom reds.

| tree | failed | passed | skipped | deselected | time |
|---|---|---|---|---|---|
| **union** `i123d-audit` + `i123d-audit-db` | **2** | **3028** | 49 | 143 | 4:23 |
| **baseline** local `plasma` + database `plasma` | **0** | **2623** | 50 | 83 | 4:07 |

**[M]** Union: `2 failed, 3028 passed, 49 skipped, 143 deselected, 54 warnings in 263.25s`,
exit 1 ([`evidence/union-unit-suite.tail.log`](i123d-audit/evidence/union-unit-suite.tail.log)).
Baseline: `2623 passed, 50 skipped, 83 deselected, 54 warnings in 247.44s`, exit 0
([`evidence/sweep/branch-plasma.suite.log`](i123d-audit/evidence/sweep/branch-plasma.suite.log)).
The baseline was run through the same stepping worktree and the same command as every branch in
§7, so it is a baseline and not a number from somewhere else.

`3028 − 2623 = 405` tests the union adds. The baseline is fully green.

### 5.1 The RMG-database suites

`python -m pytest test` from each database repository root, paired with the matching RMG-Py tree.

| tree | failed | passed |
|---|---|---|
| union `i123d-audit-db` @ `d7c8595ce` + union RMG-Py | **4** | **173** |
| baseline `plasma` @ `fb3c13c60` + baseline RMG-Py | **4** | **41** |

**[M]** `4 failed, 173 passed, 284 warnings in 22.96s` and `4 failed, 41 passed in 1.04s`
([`evidence/db-union-suite.log`](i123d-audit/evidence/db-union-suite.log),
[`evidence/db-baseline-suite.log`](i123d-audit/evidence/db-baseline-suite.log)). **The same four
failures on both bases**, all in `test/test_plasma_electron_attachment.py`:
`test_trained_species_resolve_to_their_own_rate_through_their_own_node[O2|OH|O]` and
`test_o2_rate_comes_from_the_training_set_not_a_library`. Pre-existing, identical on both, untouched
by this union, outside its scope. Reproduces the third pass's 4/173 and 4/41 exactly.

The baseline database has no `test/conftest.py` — the pin `bdbbf2864` is itself one of the merges
under audit — so its suite was run with an equivalent `rmgrc` written at that worktree's root, as
passes 1 and 3 also had to. Without that it does not fail, it *errors at fixture setup* against
whatever database `settings` resolves to.

---

## 6. Every union-only failure, named and attributed

Both union failures are the two the third pass reported. Attribution is **not** the same for the
two of them, and the §7 sweep is what separates them.

**Failure A — `test_declaration_registry_is_explicit_and_closed`. NOT union-only.** **[M]** It
fails on `i119-rr-registry` alone (`1 failed, 2684 passed`), the branch that added
`PlasmaRadiativeRecombination` to `FAMILY_ELECTRON_PLACEMENT` **[R]**
`rmgpy/electron_placement.py:172` without updating its own closure test's expected dict. The union
inherits it unchanged:

```
assert {'Plasma_Electron_Attachment': (1, 0), 'Cation_R_Recombination': (1, 0),
        'PlasmaRadiativeRecombination': (1, 0), 'PlasmaElectronImpactIonization': (1, 2)}
    == {'Plasma_Electron_Attachment': (1, 0), 'Cation_R_Recombination': (1, 0),
        'PlasmaElectronImpactIonization': (1, 2)}
```

A stale expectation in a test that exists to notice exactly this kind of addition. It is a branch
defect, not an integration defect, and it belongs to `i119-rr-registry`.

**Failure B — `test_plasma_deck_does_not_declare_family[docs/i102-quarantine/input.py]`.
Union-only.** **[M]** `i115-preflight-deck` alone is green (`2631 passed`), and `i102-quarantine`
alone is green (`2641 passed`). The test comes from the first branch; the deck it fails on comes
from the second. Neither branch can see the other, so neither could have caught it:

```
plasma deck docs/i102-quarantine/input.py declares Cation_R_Recombination, a lithium-ion-battery
SEI family excluded from every plasma configuration.
```

This is the classic integration failure — two independently-correct branches whose product is
wrong — and it is the one failure in this union that only exists because of the merge. It is also
**genuine**: `docs/i102-quarantine/input.py` really does declare the quarantined SEI family, and
`i115`'s sweep really is meant to forbid that on every plasma deck. Fixing it means editing that
deck's family list, which is out of scope here.

**No third union-only failure exists.** Beyond A and B the union's failure set is empty.

---

## 7. Each contributing branch's own state, measured this pass

Ten trees, one stepping worktree, checked out detached in turn, `check-pydas`'d, built in place,
same suite command, **serially**. Script: [`branch_sweep.sh`](i123d-audit/branch_sweep.sh); full log
[`evidence/branch-sweep.full.log`](i123d-audit/evidence/branch-sweep.full.log). Every entry printed
its resolved `rmgpy`, compiled `.so` and `database.directory` before its suite ran.

| ref | tip | database paired with | build | failed | passed |
|---|---|---|---|---|---|
| `plasma` (baseline) | `a61dc1303` | `plasma` | 0 | **0** | 2623 |
| `i110-make-guard` | `fa042c08d` | `plasma` | 0 | **0** | 2635 |
| `i119-rr-registry` | `92e2d4234` | `plasma` | 0 | **1** | 2684 |
| `i115-preflight-deck` | `ceb9fd5f1` | `plasma` | 0 | **0** | 2631 |
| `i112-marcus-work-terms` | `ab3d84072` | `plasma` | 0 | **0** | 2667 |
| `i102-quarantine` | `99602b099` | `plasma` | 0 | **0** | 2641 |
| `i126-chemkin-electrons` | `1b79bee5b` | `i123-integration-db` | 0 | **2** | 2834 |
| `i134-duplicate-drops-sink` | `9b8cc38e7` | union db | 0 | **2** | 3003 |
| `i135-tdep-roundtrip` | `948a783fd` | union db | 0 | **2** | 2853 |
| `i145-units-laundering` | `e5e42ae41` | union db | 0 | **2** | 3028 |

**[M]** Two branches are red on their own tips: `i119-rr-registry` (failure A, its own) and, from
`i126-chemkin-electrons` onward, A and B together — but `i126`, `i134`, `i135` and `i145` are tips
*on the union's own first-parent chain*, so "red on its own branch" there means "the union was
already red at that point", not "an independent branch ships red". **The only genuinely independent
branch that ships red is `i119-rr-registry`.** The third pass's finding that one branch was red on
its own branch reproduces, and it is the same branch.

The build is incremental (`make build`, no `make clean`) exactly as the third pass's sweep was, so
the counts are comparable with that pass's. `git checkout` refreshes the mtime of every file it
changes, which is what `setup.py`'s dependency check reads. **Limitation, stated:** an incremental
build cannot rule out a stale `.so` from a file that changed without its mtime moving; a full
`make clean` per branch would, at roughly triple the wall clock.

---

## 8. The lithium charge network, end to end, stage by stage

Two independent drives of the same network: a probe that calls the loader, the resolver and the
reactor directly, and the deck that drives the identical chemistry through `rmg.py` and an input
file. Both are needed — the probe proves each stage in isolation, the deck proves the input-file
path that carries them.

### 8.1 Stages 1–4, from the loader, the resolver and the reactor

**[M]** `python docs/i123-integration/probe_lithium_charge_network.py` → exit 0, every stage PASS
([`evidence/probe-lithium-network.log`](i123d-audit/evidence/probe-lithium-network.log)).

| stage | what it establishes | measured |
|---|---|---|
| 1 LOAD | both libraries load from the union database; one entry each | `[Li] => [Lip]` `VoronovEIArrhenius` electrons **+1**; `[Lip] => [Li]` `BadnellRRArrhenius` electrons **−1**; source produces exactly what the sink consumes |
| 2 BALANCE | each canonical reaction closes its own charge gap | both `is_balanced` → True |
| 3 PLACEMENT | the resolver restores the electron participants from the declaration | ionisation → `[Li] + e => [Lip] + e + e`, counts (1,2); recombination → `[Lip] + e => [Li]`, counts (1,0); both second order; canonical reactions unmutated |
| 4 REACTOR | one `PlasmaReactor` holds source **and** sink and evaluates each at its own rate law | `kf[source] = 5.031584e+07` = Voronov; `kf[sink] = 1.468188e+05` = Badnell; `electron_index = 2`; **the cation has a non-zero loss channel** |

The registry read back at stage 3 **[M]**:
`{'Plasma_Electron_Attachment': (1, 0), 'Cation_R_Recombination': (1, 0),
'PlasmaRadiativeRecombination': (1, 0), 'PlasmaElectronImpactIonization': (1, 2)}`.

### 8.2 Stages 5 and 6 — presence in the generated model, and the file

**[M]** `python rmg.py docs/i123-integration/input.py`. The composition is neutral — no cation is
seeded — so `[Lip]` has to be *produced*:

```
At time 0.0000e+00 s, species [Lip](3) was added to model core to avoid singularity
Moved 2 reactions from edge to core
    Li(2) => [Lip](3)
    [Lip](3) => Li(2)
...
The final model core has 7 species and 2 reactions
```

**Both channels are in the generated core, not merely loadable.** The written Chemkin file
([`evidence/deck-chem.inp`](i123d-audit/evidence/deck-chem.inp)):

```
REACTIONS    KCAL/MOLE   MOLES

Li(2)+e-(1)=>[Lip](3)+e-(1)+e-(1)                   1.490e+18 -0.267    162.150
    TDEP/e-(1)/   ! VoronovEIArrhenius exported as a modified-Arrhenius fit of k(Te) over 1.16e+04-2.321e+08 K; ...

[Lip](3)+e-(1)=>Li(2)                               1.734e+14 -0.801    0.010
    TDEP/e-(1)/   ! BadnellRRArrhenius exported as a modified-Arrhenius fit of k(Te) over 10-1e+07 K; ...
END
```

Both electrons stand where the reactor puts them, and both triples are **character-for-character
identical to the third pass's** — confirming the units repair moved no generated number.

**`### RMG_EXIT=1 ###`**, by design, from the Chemkin→Cantera translation failure. §10.3.

### 8.3 Stage 7–9 — three successive round trips, rate compared at each

This is where the third pass's blocker lived, so it is driven explicitly. **[M]**
`python docs/i145-units-laundering/probe_units_pairing.py <deck-dir>` → exit 0
([`evidence/probe-units-pairing.log`](i123d-audit/evidence/probe-units-pairing.log)).

|  | header | ionisation A | k(Te) source | k(Te) sink |
|---|---|---|---|---|
| as generated | `REACTIONS KCAL/MOLE MOLES` | `1.490e+18` | `1.081914e+08` | `9.618226e+04` |
| after trip 1 | same | `1.490e+18` | `1.081914e+08` | `9.618226e+04` |
| after trip 2 | same | `1.490e+18` | `1.081914e+08` | `9.618226e+04` |
| after trip 3 | same | `1.490e+18` | `1.081914e+08` | `9.618226e+04` |

**Ratios `1.000000000e+00` at every trip, for both channels.** Before the repair each trip
multiplied by 10⁻⁶ and the third file was 10⁻¹⁸ of the first. Every file in the chain declares the
same units header — the check that forbids the wrong repair, since relabelling the header to match
wrong numbers would satisfy the rate check alone.

The probe also re-establishes, on this tree, **which half was wrong**, so the repair's direction is
a measurement and not an assumption: the reader reproduces the file's printed triple character for
character and the file carries the generator's triple character for character; the residual against
the *unrounded* generator is 0.23 % and 0.38 %, both inside the 0.47 % bound that printing `n` to
three decimals sets. And the writer's reduction is now unit-independent — the same rate law declared
in `m^3/(mol*s)` and in `cm^3/(mol*s)` reduces to the identical number, ratio `1.000000e+00`.

**A ratio that is nearly one is a finding.** The probe takes `R` from `rmgpy.constants` rather than
hardcoding it; RMG carries the 2010 CODATA `8.314472 J/(mol*K)`, and hardcoding the 2018 value
makes the hand-arithmetic check fail by `1.0000079`, which is the two constants' disagreement and
nothing about the reader.

### 8.4 `scripts/mergeModels.py`, unmodified, on two half-mechanisms

**[M]** [`probe_merge_models.py`](i123d-audit/probe_merge_models.py) → exit 0
([`evidence/probe-merge-models.log`](i123d-audit/evidence/probe-merge-models.log)). The deck output
is split textually into one reaction per half — everything above the `REACTIONS` block copied
verbatim, so no number is touched on the way in — and the shipped script is invoked **as a
subprocess, by path**, as a user would:

```
Added 1 out of 1 (100.0%) unique reactions from model #1.
Added 1 out of 1 (100.0%) unique reactions from model #2.
The merged model has 7 species and 2 reactions
```

| | k(Te = 11604.5 K) |
|---|---|
| deck `Li(2)+e-(1)=>[Lip](3)+e-(1)+e-(1)` | `1.0819136429e+08` |
| merged, same reaction | `1.0819136429e+08` (**ratio 1**) |
| deck `[Lip](3)+e-(1)=>Li(2)` | `9.6182263305e+04` |
| merged, same reaction | `9.6182263305e+04` (**ratio 1**) |

**Both channels survive and neither rate constant moves.** This is the invocation through which the
third pass's 10⁶ laundering was reachable with no plasma-specific flag; it is clean now.

**One trap found writing this probe, recorded because the next person will hit it.** A Chemkin
`TDEP/` reaction reads back as a `TwoTemperaturePlasma`, whose electron-temperature entry point is
`get_rate_coefficient_two_temp(T, Te)` **[R]** `rmgpy/kinetics/arrhenius.pyx:428`. It has **no**
`get_rate_coefficient_electron_temp` — that is the spelling the two *library* rate laws use **[R]**
`arrhenius.pyx:1008,1678`. Probing only for the library spelling falls silently through to
`get_rate_coefficient(T)`, i.e. the **gas** temperature, giving `8.6e-25` for the ionisation against
`1.08e+08` at Te. The ratios still come out right, because the mistake is applied to both sides —
which is exactly why it is dangerous: the comparison passes and the printed number is 33 orders of
magnitude from what the reactor uses. Fixed in the probe and documented in it.

---

## 9. **THE BLOCKING FAILURE** — RMG's own seed mechanism cannot restart a plasma model

Every RMG run writes a `seed/` directory and a ready-to-run `restart_from_seed.py` beside its
Chemkin output, unconditionally and without being asked **[R]** `rmgpy/rmg/main.py:1745-1753,1776`.
That is a **third serialisation** of the same model, alongside Chemkin and Cantera, and no pass of
this audit has exercised it. It is reachable for a plasma model for the first time **because of this
union** — before `i126` a plasma run could not complete its first save at all.

### 9.1 What happens

**[M]** `python rmg.py <deck>/restart_from_seed.py`, the file RMG itself generated, unedited:

```
Adding seed mechanism restart to model core...
Added 2 new core reactions
    Li(2) => [Lip](3)
    [Lip](3) => Li(2)

Adding reaction library PlasmaElectronImpactIonization to model edge...
Warning: Marked reaction Li(2) => [Lip](3) as duplicate of Li(2) => [Lip](3) for saving to Chemkin file.
    Created 1 new edge reactions
        Li(2) => [Lip](3)

Adding reaction library PlasmaRadiativeRecombination to model edge...
This library reaction was not new: [Lip](3) => Li(2)
...
    The model core has 7 species and 3 reactions
```

**Three reactions where the original run had two.** The ionisation channel is in the core **twice**;
the recombination once. Then, at the first Chemkin save:

```
rmgpy.exceptions.MechanismWriterError: Reaction Li(2) => [Lip](3) has kinetics VoronovEIArrhenius
whose rate is proportional to the electron density, but the exported equation
"Li(2)=>[Lip](3)+e-(1)" has no electron among its reactants, so a solver would evaluate it at the
wrong reaction order.
### RESTART_EXIT=1 ###
```

[`evidence/restart-plasma-with-libraries.stdout.log`](i123d-audit/evidence/restart-plasma-with-libraries.stdout.log),
[`.stderr.tail.log`](i123d-audit/evidence/restart-plasma-with-libraries.stderr.tail.log).

### 9.2 Isolated: it is the seed, not the duplication

Re-run with `reactionLibraries=[]` so the seed is the only source of chemistry **[M]**: core is 7
species and **2** reactions — the right model — and it **still dies at the same save, same
exception, exit 1**
([`evidence/restart-plasma-seed-only.stdout.log`](i123d-audit/evidence/restart-plasma-seed-only.stdout.log)).
So the write refusal is caused by the seed alone; the doubled channel is a second, separate
consequence of the same root cause.

### 9.3 The root cause, derived and then measured

**The electron-placement declaration is keyed on the reaction's owner label, and a seed mechanism
does not preserve the owner.**

- **[R]** `rmgpy/rmg/main.py:1881` — `kinetics_library.name = "seed"`, with the comment *"Rename for
  the output directory, as these names should not be dynamic"*. The originating library's label is
  discarded here; it survives only as free text in `longDesc`, which nothing parses.
- **[R]** `rmgpy/rmg/main.py:646,657-658` — on restart the directory is copied to `restart/` and
  appended to `seed_mechanisms` under the label `"restart"`.
- **[R]** `rmgpy/electron_placement.py:172` — `FAMILY_ELECTRON_PLACEMENT` is a dict keyed by owner
  label. Neither `"seed"` nor `"restart"` is a key.
- **[R]** `rmgpy/electron_balance.py:135` — with no declaration the counts fall back to the
  net-derived rule: `electrons < 0 → (−electrons, 0)`, `electrons > 0 → (0, electrons)`.

**[M]** [`probe_seed_roundtrip.py`](i123d-audit/probe_seed_roundtrip.py) → exit 0
([`evidence/probe-seed-roundtrip.log`](i123d-audit/evidence/probe-seed-roundtrip.log)) measures each
step of that on the objects themselves:

| | canonical | via the seed | consequence |
|---|---|---|---|
| electron-impact ionisation | placement **(1, 2)** | placement **(0, 1)** | **declaration lost** |
| radiative recombination | placement **(1, 0)** | placement **(1, 0)** | preserved *by luck* |

The seed is otherwise **lossless**: it keeps the original `VoronovEIArrhenius` and
`BadnellRRArrhenius` objects, not a Chemkin reduction, and keeps the scalar electron counts `+1` and
`−1`. Only the owner is gone.

**That table explains the asymmetry in §9.1 exactly.** For the recombination (`electrons = −1`) the
fallback gives `(1, 0)`, which is what the declaration says anyway — so the seed copy and the
library copy compare as the same reaction and RMG says *"This library reaction was not new"*. For
the ionisation (`electrons = +1`) the fallback gives `(0, 1)` where the declaration says `(1, 2)` —
the **consumed** electron vanishes from the reactant side, which is precisely the incident-order
information a net scalar cannot carry — so the two copies compare as different reactions and RMG
adds the channel a second time and marks it `DUPLICATE`. Measured directly **[M]**:
`is_isomorphic(either_direction=True)` is `False` for the ionisation pair and `True` for the
recombination pair.

### 9.4 The negative control: an ordinary mechanism restarts perfectly

**[M]** `examples/rmg/minimal`, same RMG-Py tree, same database, run then restarted from its own
generated `restart_from_seed.py`
([`evidence/restart-neutral-control.log`](i123d-audit/evidence/restart-neutral-control.log)):

| | final core | duplicate warnings | exit |
|---|---|---|---|
| original run | 26 species, **66 reactions** | 0 | 0 |
| restart from its own seed | 26 species, **66 reactions** | 0 | **0** |
| plasma deck original | 7 species, **2 reactions** | 0 | 1 (Cantera, §10.3) |
| plasma deck restart | 7 species, **3 reactions** | 1 | **1 (write refused)** |

The restart machinery is not broken. **It is broken for the electron representation this union
delivers.**

### 9.5 Why this blocks, and the argument against, stated fairly

**No wrong number ever reaches a file or a solver.** **[R]** `rmgpy/rmg/main.py:967` —
`save_everything()` runs immediately after the initial enlarge and **before** the model-generation
loop that simulates. The doubled model is therefore refused at the writer before any solver sees it,
deterministically, on every path. The failure is loud, exit 1, with an actionable message.

**Why I still call it blocking:**

1. **It is the same shape as the campaign's own first blocker.** Pass 1 blocked on *"the Chemkin
   writer refused every plasma rate law, so a generated mechanism could not be saved."* This is that
   sentence with "generated" replaced by "restarted". Loudness did not make that one acceptable.
2. **There is no workaround for the use case.** The two known non-blocking defects each have one —
   turn the RMS writer off in the deck (§10.1), split the channels across two libraries (§10.2). A
   user who wants to restart a long plasma run has none short of hand-editing the seed's library
   label, and RMG offers no way to decline writing a seed.
3. **The artifact is produced silently.** RMG prints `Making seed mechanism...` and says nothing
   about the seed being unusable. The loudness arrives in a *different run*, possibly days later.
   This campaign has already ruled, in `36a43f54c` *"Fail the run when an end-of-run export step
   produces nothing"*, that quietly emitting a broken artifact is not acceptable.
4. **It is a third instance of one structural weakness**, not a new one: the placement declaration is
   keyed one-placement-per-owner-label and does not survive serialisation. Pass 1's blocker, the
   strict-`xfail` landmine in §10.2, and this are the same obstruction reached from three
   directions. Whoever fixes it should fix it once.

**The argument against, which the manager may reasonably take:** restart-from-seed is not among the
paths this brief enumerates as the deliverable, the failure emits no wrong numbers, and merging this
union breaks nothing that works on `plasma` today — on `plasma` a plasma deck cannot run at all. On
that reading this is a **fifth known non-blocking defect** and the verdict is READY. For that
reclassification to be safe it needs, at minimum: a ticket owning the owner-preserving fix, and a
loud warning at seed-write time when the model contains any reaction whose owner carries a placement
declaration — so that the failure is announced by the run that causes it rather than by the run that
suffers it.

I have not made that change: this pass fixes nothing, by instruction.

---

## 10. The four known defects — confirmed, not re-litigated, not fixed

**[M]** `python docs/i123c-audit/probe_known_defects.py` → exit 0
([`evidence/probe-known-defects.log`](i123d-audit/evidence/probe-known-defects.log)).

### 10.1 The RMS YAML writer still refuses these rate laws — still **loudly**

**[R]** `rmgpy/yaml_rms.py:252`, and the preceding line is `else:` — the dispatch's terminal branch
is a `raise`, not a fall-through default, so an unknown rate law cannot be silently mis-emitted.
Driven on both objects **[M]**:

```
VoronovEIArrhenius  -> ValueError: Object of type <class '...VoronovEIArrhenius'> does not have a
                       defined conversion to ReactionMechanismSimulator format
BadnellRRArrhenius  -> ValueError: (same shape)
```

The property the decision rests on holds. The shipped deck still carries `generateRMSYAML=False` as
its documented workaround.

### 10.2 One library carrying both channels still cannot be loaded

**[M]** `check_for_duplicates` refuses:
`DatabaseError: Unexpected duplicate reaction [Li] => [Lip] in kinetics library PlasmaBothChannels.
Reaction index 2 matches index 1.` The two entries carry `family=None` at that point — the owner is
attached later, in `get_library_reactions`, which is why the placement declaration is invisible to
the duplicate check.

**The pin is still expected-*failing*, not silently passing** **[M]**
([`evidence/xfail-pin.log`](i123d-audit/evidence/xfail-pin.log)):

```
XFAIL test/rmgpy/i134DuplicateElectronsTest.py::TestTheLibraryLoadLandmineIsNotReachedFromHere::
      test_one_library_carrying_both_channels_can_be_loaded
======================== 179 passed, 1 xfailed in 1.78s ========================
```

`@pytest.mark.xfail(strict=True, ...)` **[R]** `test/rmgpy/i134DuplicateElectronsTest.py:582`, so an
XPASS would fail the suite. The shipped deck still keeps the two channels in two libraries.

### 10.3 The Chemkin→Cantera translation failure now surfaces in the exit status

**[M]** The deck run ends with:

```
EXPORT FAILED. 1 of RMG's output steps did not produce their files:
  Chemkin-to-Cantera translation (cantera_from_ck/)
      InputError: Error while reading reaction in chem.inp starting on line 66:
### RMG_EXIT=1 ###
```

It no longer exits 0. The underlying cause is the ecosystem converter, §10.5.

### 10.4 The lithium cation's enthalpy — **measured on this union**

Not inherited from either the old description or the fix. **[D]** `LithiumPrimaryThermo`, entry
`[Lip]`, index 65 **[M]**:

```
H298      = 532.936 kJ/mol
E0        = 526.738 kJ/mol
shortDesc = []        <- empty: no provenance recorded
reference dfH(Li+, g) = dfH(Li, g) + IE(Li) = 159.3 + 520.2 = 679.5 kJ/mol
gap = 146.6 kJ/mol (1.52 eV) low
```

**[I]** The reference is the standard atomisation enthalpy plus the first ionisation energy of
lithium, quoted rather than re-derived from a primary source here. `i129-lithium-cation-enthalpy` @
`75e0cb6f0` is **not** an ancestor of this union (§2.2), so the gap is exactly as the third pass
measured it. The lithium libraries present are `LithiumPrimaryThermo`, `LithiumAdditionalThermo`,
`PlasmaCationThermo`.

**Consequence, bounded:** both channels in the deck are irreversible, so no `Keq` is formed from
this enthalpy and no rate in §8 depends on it. It would matter to any reversible cation chemistry
and to any equilibrium composition.

### 10.5 The ecosystem Chemkin converter — the smaller claim, said small

Cantera's `ck2yaml` does not implement the `TDEP/` auxiliary keyword RMG's writer emits **[M]**:

```
Error while reading reaction in chem.inp starting on line 66:
"Li(2)+e-(1)=>[Lip](3)+e-(1)+e-(1)   1.490e+18 -0.267  162.150
     TDEP/e-(1)/ ..."
could not convert string to float: 'e-(1)'
```

from `/home/alon/Code/Cantera/build/python/cantera/ck2yaml.py:2216`. **"Round-trips" in this report
means RMG reads back what RMG wrote.** The file does not import into Cantera, and nothing here
should be read as claiming otherwise. RMG's own Cantera writer (`generateCanteraYAML2=True`) is a
separate path and did produce its files.

---

## 11. The reloaded Chemkin representation — a property, not a defect

**[M]** `python docs/i123c-audit/probe_readback_identity.py <deck-dir>` → exit 0
([`evidence/probe-readback-identity.log`](i123d-audit/evidence/probe-readback-identity.log)).

A reaction read back from Chemkin carries its electrons as **explicit species in the participant
lists with a zero scalar count**, where the canonical form carries a scalar and no electron
participant:

| | canonical | reloaded |
|---|---|---|
| type | `LibraryReaction` | `LibraryReaction` |
| kinetics | `VoronovEIArrhenius` / `BadnellRRArrhenius` | `TwoTemperaturePlasma` (both) |
| `electrons` scalar | `+1` / `−1` | `0` / `0` |
| placement counts | `(1, 2)` / `(1, 0)` | `(0, 0)` / `(0, 0)` |
| electron participants | none | 1 in / 2 out, and 1 in / 0 out |
| owner | the two library labels | `Unclassified` |

**Q1 — is a reloaded reaction the same reaction as its canonical original?** **No**, in either
direction, for either channel **[M]**. The round trip recovers the **chemistry** in the expanded
electron representation, not the canonical **object**. A reloaded mechanism and the model that
produced it are not interchangeable under RMG's identity predicate.

**Q2 — do the two reloaded channels collapse into each other** (the I-134 defect seen from the far
side of the file)? **No** **[M]**. In the expanded form the electrons are ordinary participants, so
the two channels differ in their participant lists alone: the expansion is **self-distinguishing**
and needs no declaration.

**Q3 — can a reloaded plasma mechanism be fed back into a model?** **No, and loudly** **[M]**. Both
routes refuse: `Reaction.is_balanced` → `False` for both (the electron is counted as a conserved
element on both sides while the scalar is zero, so the tallies differ and the scalar cannot
compensate), and `CoreEdgeReactionModel.make_new_reaction` raises `Could not retrieve the
family/library: Unclassified`.

**Reported as a property**, per the brief, because every consequence of it is a refusal rather than
a wrong number. It is worth stating plainly what the property is: **Chemkin is an export format for
this chemistry, not an interchange format.** §9 is what happens when a *different* serialisation —
one that keeps the canonical form and looks like an interchange format — loses one field.

### 11.1 The three consumers of the repaired identity predicate, re-measured

**[M]** `python docs/i123c-audit/probe_identity_consumers.py <deck-dir>` → exit 0
([`evidence/probe-identity-consumers.log`](i123d-audit/evidence/probe-identity-consumers.log)). All
green on this tree: `ReactionModel.merge` keeps both channels and both controls still collapse what
should collapse; a `PDepReaction` is born with `electrons = 0` and placement `(0,0)` and does not
match either charged channel; and the products-only generation shortcut is reached on a plasma path,
carrying charged reactions, with the doubted negation never taken.

---

## 12. The negative control

`examples/rmg/minimal`, run three ways to separate the two variables **[M]**
([`evidence/negative-control.log`](i123d-audit/evidence/negative-control.log)):

| tree | database | exit | `chem.inp` sha256 |
|---|---|---|---|
| union `i123d-audit` | baseline `plasma` | **0** | `b318e4c74f03c5fcfbf2d1c13ce28d1076cf17e4bd5619855c2bf6cb7c317fb8` |
| baseline `plasma` | baseline `plasma` | **0** | `b318e4c74f03c5fcfbf2d1c13ce28d1076cf17e4bd5619855c2bf6cb7c317fb8` |
| union `i123d-audit` | union database | **0** | `b318e4c74f03c5fcfbf2d1c13ce28d1076cf17e4bd5619855c2bf6cb7c317fb8` |

`chem_annotated.inp` → `2850089d64e3d74423d22bc26ad2b1969741c93256ebd00beacb34a99cdb54e5` and
`species_dictionary.txt` → `5114ffb2f4b357170486a77d5b0ad86808eac4b25d434665e9cbeaf76bd85e8f`, all
three identical across all three runs, and `diff -r -q` over the whole `chemkin/` directory of the
first two exits **0** with no output.

The third row is the one worth having: holding the tree fixed and swapping the database changes
nothing either, so the union database's extra libraries and argon thermo are inert for ordinary
chemistry. **Same bytes, same exit status, both variables isolated.**

The independent byte control inside the units probe agrees: an ordinary non-plasma mechanism written
three times is byte-identical each time, `4021 bytes`, digest
`ba8e2c58cea509e93069d432f47b1d3808cf0d59114495fa635f896a5ac7c65c`.

---

## 13. Verdict

# NOT READY

**The blocking failure:** RMG's seed mechanism does not preserve a reaction's owner
(`rmgpy/rmg/main.py:1881`), the electron-placement declaration is keyed on exactly that owner
(`rmgpy/electron_placement.py:172`), and so **every plasma reaction that passes through a seed
silently loses its declaration**. The consequences, both measured on this union: the generated
`restart_from_seed.py` builds a core carrying the ionisation channel **twice**, and the run then
**dies at its first Chemkin save with `MechanismWriterError`, exit 1**. With the reaction libraries
removed the model is right and the save still fails, so the seed alone is sufficient to cause it. An
ordinary non-plasma mechanism restarts from its own seed to a byte-identical model, exit 0. §9.

**What is genuinely green, and is not in doubt:**

- All three previous blockers are repaired on this base and stay repaired under exercise: the writer
  emits both plasma rate laws (§8.2), the sink is not discarded (§8.1, §8.2), and three successive
  round trips move no rate constant by a digit (§8.3).
- The lithium charge network runs end to end, stage by stage, twice over — through the API and
  through an input file (§8).
- `scripts/mergeModels.py`, unmodified, keeps both channels with both rate constants exact (§8.4).
- The negative control is byte-identical to the shared branch with both variables isolated (§12).
- The baseline shared branch is fully green, 2623 passed / 0 failed, same command (§5).

**Also true, and not blocking:**

- The branch is **not merge-clean**: four `rmgrc` pin commits must come out first (§3).
- Two unit failures: one inherited from `i119-rr-registry`, one genuinely union-only from the
  `i115` × `i102` pairing (§6).
- Four known defects confirmed as described, with the lithium cation enthalpy measured at
  **146.6 kJ/mol low** on this union (§10).

**If the manager reclassifies §9 as a fifth known non-blocking defect** — a defensible reading,
argued fairly in §9.5 — then every other Verifier item is satisfied and the answer becomes READY
subject to the `rmgrc` reverts. That decision is not mine and not the manager's alone; §9.5 names
the two conditions that would make it safe. **It was put, and it was decided — §13.1.**

### 13.1 Both open questions were put to the operator, and both were decided

Recorded here so the next reader finds a settled question rather than an auditor's judgement call
left hanging. Put after the report was written and committed, with §9.5's counter-argument stated
as its own option rather than as a footnote to the recommendation:

| question | decision |
|---|---|
| Does §9 block, or is it a fifth known non-blocking defect? | **Stands as blocking.** NOT READY is upheld. |
| Who removes the four `rmgrc` pins, and when? | **Left in place; whoever performs the merge reverts them**, as the first step. |

Of the four reasons in §9.5, two are the load-bearing pair behind the ruling: **the defect has no
workaround** — the RMS writer can be switched off and the two channels can be split across two
libraries, but a seed cannot be declined — and **RMG emits the broken artifact silently**, so the
failure surfaces in a different run, possibly days later, rather than in the run that caused it.

The pins stay because this branch needs them to load its own database: removing them would make
every command in §15 read the wrong tree, and the report would stop being reproducible from its own
tip. §3 carries all four SHAs — `ae52de38a`, `9874fc706`, `34de41c76`, `cd6cde425` — and the line
they must be restored to. **That is a merge-time obligation on the merger, not an outstanding task
of this audit.** It is also the fourth pass running to inherit the previous pass's pin by not
reading it (§0.2) — a recurring defect in this lineage's configuration rather than four independent
slips. Making `rmgrc` un-inheritable was considered as part of the same decision and deliberately
left out of scope here; it wants its own ticket.

---

## 14. What this pass could not reach

1. **Whether the doubled ionisation channel can ever be simulated.** Argued from `main.py:967` —
   `save_everything()` precedes the model-generation loop unconditionally — and observed on two
   runs, but not proven by exhausting every code path into the solver. If a configuration exists
   that simulates before its first save, §9 becomes a silent-wrong-number defect rather than a loud
   one.
2. **The other two seed artifacts.** `seed_edge/` and `filters/` were not exercised; only the core
   seed was. A plasma deck with a non-empty edge may fail differently or earlier.
3. **Any plasma chemistry other than lithium.** Argon and the electron-attachment network were not
   driven end to end. The four pre-existing `test_plasma_electron_attachment` failures (§5.1) sit in
   exactly that area and were not investigated — they are identical on the baseline, so they are not
   this union's, but "not this union's" is not "understood".
4. **Whether `i119-rr-registry`'s red (failure A) has any consequence beyond the stale expectation.**
   I read the assertion and stopped there, as instructed not to fix.
5. **The incremental-build caveat in §7.** No branch in the sweep got a `make clean`; a stale `.so`
   from a file whose mtime did not move would not have been caught.
6. **RMS end to end.** Verified at the writer and by reading the dispatch's terminal branch, not by
   running a job with `generateRMSYAML=True` to see the failure arrive at save.
7. **The enthalpy reference.** Quoted (159.3 + 520.2 kJ/mol), not re-derived from a primary source.

---

## 15. Reproducing this

```bash
conda activate rmg_env
cd /home/alon/Code/RMG-Py-i123d-audit
export PYTHONPATH=/home/alon/Code/RMG-Py-i123d-audit:$PYTHONPATH
make build                                   # never bare `make`

python -m pytest -m "not functional and not database"          # union suite
bash docs/i123d-audit/branch_sweep.sh                          # baseline + 9 branches, serial

DECK=<scratch>/deck && mkdir -p $DECK && cp docs/i123-integration/input.py $DECK/
python rmg.py $DECK/input.py                                   # exits 1 by design, §10.3

python docs/i123-integration/probe_lithium_charge_network.py   # stages 1-4
python docs/i145-units-laundering/probe_units_pairing.py $DECK # three round trips
python docs/i123d-audit/probe_merge_models.py       $DECK      # mergeModels
python docs/i123c-audit/probe_readback_identity.py  $DECK      # reloaded representation
python docs/i123c-audit/probe_identity_consumers.py $DECK      # identity consumers
python docs/i123c-audit/probe_known_defects.py                 # the four known defects
python docs/i123d-audit/probe_seed_roundtrip.py     $DECK      # the blocking failure, §9

cp -r $DECK/seed <scratch>/restart/ && cp $DECK/restart_from_seed.py <scratch>/restart/
python rmg.py <scratch>/restart/restart_from_seed.py           # exits 1: THE BLOCKER
```

### 15.1 The scratch worktrees this pass created

All four are **detached**, carry no branch, and hold nobody's work but this audit's. They are
**left in place**, not deleted — §0.2 is what happens when a later pass silently inherits a pointer
to a worktree it does not own, and a named, detached, documented worktree is the opposite of that.
Named here so nobody mistakes them for someone's live work:

| path | ref | what it was for |
|---|---|---|
| `/home/alon/Code/RMG-Py-i123d-step` | detached; ends at `e5e42ae41` | stepped through the ten refs of §7 |
| `/home/alon/Code/RMG-Py-i123d-baseline` | detached at `plasma` `a61dc1303` | the §12 negative control's baseline arm |
| `/home/alon/Code/RMG-database-i123d-baseline` | detached at database `plasma` `fb3c13c60` | baseline database for §5.1 and §7 |
| `/home/alon/Code/RMG-database-i123d-step1` | detached at `i123-integration-db` `c7bd96292` | the database `i126-chemkin-electrons` was built against |

Each carries an **untracked** `rmgrc` written by the sweep or by hand; none is committed anywhere.
To remove them once this pass is closed:

```bash
git -C /home/alon/Code/RMG-Py-i123d-audit worktree remove --force /home/alon/Code/RMG-Py-i123d-step
git -C /home/alon/Code/RMG-Py-i123d-audit worktree remove --force /home/alon/Code/RMG-Py-i123d-baseline
git -C /home/alon/Code/RMG-database-i123d-audit worktree remove --force /home/alon/Code/RMG-database-i123d-baseline
git -C /home/alon/Code/RMG-database-i123d-audit worktree remove --force /home/alon/Code/RMG-database-i123d-step1
```

No shared checkout was written to. Nothing was merged, pushed, rebased, amended, force-pushed or
deleted, in either repository.
