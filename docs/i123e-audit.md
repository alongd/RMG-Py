# I-123 fifth audit (i123e) — plasma/electron-chemistry union

Verdict: **PENDING** (draft; filled as measurement completes).

Evidence labels: **[R]** code (file:line) · **[D]** database · **[M]** measured (command + output) ·
**[I]** inference with basis.

Worktrees:
- RMG-Py: `/home/alon/Code/RMG-Py-i123e-audit` @ branch `i123e-audit`
- RMG-database: `/home/alon/Code/RMG-database-i123e-audit` @ branch `i123e-audit-db`

Environment confirmed [M]:
- `python -c "import rmgpy; print(rmgpy.__file__)"` → `/home/alon/Code/RMG-Py-i123e-audit/rmgpy/__init__.py`
- resolved `database.directory` → `/home/alon/Code/RMG-database-i123e-audit/input` (exists). Printed
  at the head of every probe below.

---

## 1 + 16. Merge inventory (shown from git)

### RMG-Py `i123e-audit` — 46 commits over local `plasma` (merge-base `a61dc1303`) [M]
`git log --oneline plasma..i123e-audit` (46 commits). Notable contributing merges/commits:
- `cc508be3e` rmgrc repoint (this pass, first commit) — pin #6
- `278abfe32 / 69cf533f7 / 308877556` I-148 seed-placement repair + warning + round-trip pin
- `ce93e6c12` rmgrc pin #5 (I-148)
- `cd6cde425` pin #4, `34de41c76` pin #3, `9874fc706` pin #2, `ae52de38a` pin #1
- `1b7bafc45` merge I-135 (Chemkin TDEP reader repair): `36a43f54c`, `89c896a19`
- `7758c27d2` merge I-134 (reaction-identity by electron placement): `b8204f68b`
- `1b79bee5b` export path reads electron-placement declaration
- `f087cdf04` merge i102-quarantine; `5e1686bcb` merge i112-marcus-work-terms;
  `9f0eae40b` merge i115-preflight-deck; `c1f23e93e` merge i119-rr-registry;
  `582fd6ab9` merge i110-make-guard
- electron-representation commits (`0af07c2a2`, `98b0e70b4`, `aee9e510e`, `04bbb2881`, ...)

### RMG-database `i123e-audit-db` — 21 commits over local `plasma` (merge-base `fb3c13c60`) [M]
`git log --oneline plasma..i123e-audit-db` (21 commits). Notable:
- `2d7123b08` argon cation thermochemistry (sourced)
- `c7bd96292` merge i119-recombination (`9325c2d27` radiative recombination library owner)
- `1fb224371` merge i114-ionisation-declaration (`6862c21c5` electron-impact ionisation library owner)
- `bdbbf2864` test: pin database.directory for the suites in this repository
- `ac944618e` merge i104-alkali-plasma; `79320d363` merge i103-electrochem-provenance;
  `3b61056c1` merge i111-sei-reclassification (`03b96a974` reclassify, `7d4f47d3a` quarantine Marcus data)

## 3. Local-vs-remote base gap (verdict-neutral) [M]

- RMG-Py: local `plasma` = `a61dc1303`, `origin/plasma` = `6d3c03823` → local **10 ahead, 0 behind**.
  The 10 local-only commits are the walls A–D / Marcus-dGrxn / I-098 reverse-reconstruction work,
  including `fa3b58f05` which pins rmgrc to `../RMG-database-plasma/input`. The union's merge-base
  with local `plasma` **is** local `plasma`'s tip (`a61dc1303`), so testing against local `plasma`
  is testing against the true base. Remote is stale; using it would test against a base the union was
  never built on. Gap is verdict-neutral: it is entirely "local has more," none on the union's side.
- RMG-database: local `plasma` = `origin/plasma` = `fb3c13c60` → **0 ahead, 0 behind**. No gap.

## 4. Build [M]

- `python utilities.py check-pydas` then `make build` → **exit 0** [M] (`BUILD_EXIT=0`; incremental
  in-place Cython build of the full ext_modules list, no errors).
- Bare `make` still refused [M]: exits 2 with "Refusing to modify the shared RMG environment. / Use
  `make build` for an in-place worktree build." (Makefile:72 guard). No pip invoked.

## 13. The two configuration changes — safe merge order [M]

Untracking sibling branch: **`i149-rmgrc-untracked`** (tip `cdb706509`) = local `plasma` + 2 commits:
`fccda5fc0` (git-rm `rmgrc`, add `rmgrc.template`, `.gitignore /rmgrc`, add
`Settings.require_database_directory()` so a missing/wrong `database.directory` fails loud) and
`cdb706509` (CLAUDE.md doc). `git ls-tree i149-rmgrc-untracked -- rmgrc` is empty [M].

- `git --version` = 2.43.0. `git merge-tree --write-tree --messages i149-rmgrc-untracked i123e-audit`
  → exit 1, **CONFLICT (modify/delete): rmgrc deleted in i149-rmgrc-untracked and modified in
  i123e-audit** (stage-1 base blob `90292bed` == plasma's rmgrc; stage-3 = union `33f2cec8`) [M].
- `git merge-tree ... i149-rmgrc-untracked plasma` → exit 0, clean (i149 descends plasma) [M].
- **31 local branches still track `rmgrc`** ("~30" ✓): 17 at the plasma baseline blob `90292bed`,
  14 carrying their own DB pin (incl. i123e-audit) [M].
- Reverting the union's **six** rmgrc pins restores `rmgrc` to blob `90292bed`, **byte-identical to
  plasma** (`ae52de38a^ == plasma`; only the six pins touch the file) [M].

**Safe order (recommended): land `i149-rmgrc-untracked` on `plasma` FIRST, then merge the union with
its six DB pins reverted.**
- (a) untracking-first, then union with pins reverted → union's rmgrc == base `90292bed`, delete
  applies cleanly, **no conflict** (demonstrated via baseline tracking branches i093/i143, merge-tree
  exit 0) [M]. Union with pins *intact* → modify/delete conflict, must resolve toward deletion [I].
- (b) reverse (union into plasma first with pins reverted, then untracking) → also clean [M/I].
- (c) reverting the six pins is **sufficient in either order**; the union need not additionally
  `git rm rmgrc` (i149 already does) [M/I].
- (d) a worker on one of the 14 modified branches gets a modify/delete conflict (resolve toward
  deletion → gitignored `rmgrc` from template); the 17 baseline branches merge clean. Same in both
  orders, since untracking lands on plasma regardless [M].

Reason: untracking-first makes every later branch merge against a plasma that already refuses silent
DB defaults, and since the union must revert its pins anyway (making rmgrc == base), that merge is a
clean deletion rather than a modify/delete conflict.

## 5. Full unit suites (serial, same command, [M])

Command (identical on both trees), run serially, no xdist:
`python -m pytest -m "not functional and not database" -o addopts="--keep-duplicates --ignore=test/regression" -q -p no:cacheprovider`

| Tree | rmgpy resolves to | database | result |
|---|---|---|---|
| Union `i123e-audit` @ cc508be3e | `/home/alon/Code/RMG-Py-i123e-audit/rmgpy` | `/home/alon/Code/RMG-database-i123e-audit/input` | **2 failed, 3045 passed, 49 skipped, 143 deselected** (278s), exit 1 |
| Base `plasma` @ a61dc1303 (own worktree, own build) | `/home/alon/Code/RMG-Py-i123e-baseline/rmgpy` | `/home/alon/Code/RMG-database-plasma/input` (branch plasma) | **0 failed, 2623 passed, 50 skipped, 83 deselected** (156s), exit 0 |

No concurrent pytest ran during either (checked `ps` before each; runs were serial). Base DB resolves to
`RMG-database-plasma/input` when queried directly (a transient empty print in the launcher's inline
`python -c` was a shell-quoting artifact, not a resolution failure — verdict-neutral).

## 6. Union-only failures (named) [M]

Base has **0** failures. The union has exactly **2**, both union-only, both the documented known
non-blockers (item 12.v), neither a new blocker:

1. `test/rmgpy/electronPlacementTest.py::TestElectronPlacementResolver::test_declaration_registry_is_explicit_and_closed`
   — the closed-registry pin. The test file is **present on base and passes there** (base
   `FAMILY_ELECTRON_PLACEMENT` = 2 entries, old `('reactants', 1)` format). The union grew the registry
   to 4 entries (added `PlasmaRadiativeRecombination`, `PlasmaElectronImpactIonization` via the i114/i119
   chain) and reformatted to `(reactant, product)` tuples, so the closed-set pin no longer matches. A
   test-not-updated failure; do-not-fix per brief.
2. `test/rmgpy/preflightDeckFamilyExclusionTest.py::PlasmaDeckFamilyExclusionTest::test_plasma_deck_does_not_declare_family[docs/i102-quarantine/input.py]`
   — **absent on base** (union-introduced by i115). Red because the i102-quarantine *demonstration* deck
   deliberately declares `Cation_R_Recombination` to exercise the quarantine hard-fail; the exclusion
   sweep globs all decks and meets it. Do-not-fix per brief.

## 8. The newest repair (I-148) — four questions [M/R]

Probe: `docs/i123e-audit/probe_i148_repair.py` (object-level; head prints resolved rmgpy + DB).

- **8.1 Monotonicity — HOLDS [M].** `get_placement_owner` (`rmgpy/electron_balance.py:135-141`) reads
  `family` first, `library` provenance only when `family` names no declaration. Measured across an
  ordinary H_Abstraction (family-generated), a plasma family reaction, a fresh library reaction, a
  seed-renamed reaction with preserved `library`, and a seed-renamed reaction with no provenance: every
  reaction that resolved under the pre-I-148 (family-only) rule resolves **identically**; the **only**
  newly-resolving reaction is the seed-renamed one whose `family` names nothing and whose `library`
  provenance is preserved — exactly the I-148 case. Family-generated and library reactions are untouched.
- **8.2 Guard not weakened — INTACT [M].** `resolve_electron_placement` (`rmgpy/electron_placement.py:323-328`)
  raises `ElectronPlacementError` ("carries no family attribution; electron placement is family-declared
  and cannot be inferred from the net electron count") for a plasma-shaped reaction with no family, and
  **never** falls back to the net rule. Shown.
- **8.3 Duplicated-parsing coupling — costs a wrong WARNING, not a wrong number [R].** The write-time
  predicate `seed_placement_survives` (`rmgpy/data/kinetics/library.py:307-318`) mirrors the reader
  `KineticsLibrary.get_library_reactions` (`library.py:342-365`, provenance parse `:343-345`, container
  rename `rxn.family = self.label` `:355`) in the reader's own branch order (library-provenance →
  rate-rule-family → lost). But the predicate feeds **only** `warn_if_seed_loses_placement`
  (`rmgpy/rmg/main.py:2386`, warning `SEED LOSES ELECTRON PLACEMENT` `:2419`) — an advisory diagnostic at
  seed-write time. The actual reload uses the reader + `get_placement_owner` as the source of truth. So
  drift between the two implementations produces a **false or missing warning**, never a wrong rate
  constant or wrong placement. This is a materially weaker coupling than the I-126 two-converters-drift it
  was compared to (there both implementations produced numbers); here one produces only a log line.
- **8.4 Latent seed→core defect — see item 9 seed restarts.** The net-derived fallback lives at
  `rmgpy/electron_balance.py:292-294` (and the disagree branch `:264-274`, whose docstring states no
  shipped path reaches it). The claim is that a family-less seed reaction with a nonzero scalar electron
  count and no usable provenance would demote and drop electrons. Whether the shipped seed path reaches
  it is measured by the two whole-run seed restarts in item 9 (if the restarted core carries the
  ionisation `(1,2)` intact and single-copy and saves byte-identically, the shipped path does not reach
  it). Not fixed either way.
