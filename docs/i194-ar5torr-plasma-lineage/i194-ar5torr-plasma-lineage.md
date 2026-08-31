# I-194 — the argon deck re-run on a runtime cut from the source of truth

Every argon run this campaign reported ran on `/home/alon/Code/RMG-Py-argonrun`, a runtime
**28 commits behind** `plasma`. This ticket rebuilt the runtime from source at `plasma`'s tip and
re-ran the identical `PlasmaArgon` deck to get a trustworthy measurement. **No RMG source was
changed.** The prior run's record is
`/home/alon/Code/RMG-database-i193-plasmaargon/docs/i193-plasmaargon.md`.

## Lineage

```
$ git log --oneline HEAD..plasma | wc -l
0
```

This worktree (`i194-plasma-runtime`, tip `3b479a638`) *is* `plasma`'s tip. The 28-commit gap the
prior runtime carried is closed.

## Runtime

Built from source, **not** copied `.so`:

```
$ python utilities.py check-pydas          # writes rmgpy/solver/settings.pxi (DEF DASPK = 1)
$ python setup.py build_ext --inplace -j 8
REAL_BUILD_EXIT=0
```

(The contract's bare `setup.py build_ext` fails with `'settings.pxi' not found`; `check-pydas`
generates that build file first, exactly as `make build` does. `settings.pxi` is generated, not
committed.) The built `rmgpy/reaction.cpython-39...so` has a distinct inode and a fresh 08:30
mtime, different in size from argonrun's copy — a from-source build, not a copy.

Three fixes named in the contract are all present in this runtime:
`1458a4592` (stop injecting Ar/He/Ne/N2), `ebb95ad96` (skip Chemkin→Cantera for plasma),
`ddb91501e` (opposite-charge species no longer collapse).

## The run

Deck: `docs/i194-ar5torr-plasma-lineage/input.py`, copied verbatim from
`RMG-Py-argonrun@d857d0b16:docs/i186-ar5torr-firstlight/input.py` (loads `PlasmaArgon`, no
`generatedSpeciesConstraints` directive). `rmgrc` was **absent** before this ticket (recorded:
no prior value); set to `/home/alon/Code/RMG-database-i193-plasmaargon/input`.

```
$ python rmg.py input.py > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)
PYTHON_EXIT=1
```

## The six questions

### 1. Are He, Ne, N2 still injected into the core? — **NO (was YES)**

Core is **3 species / 1 reaction**: `e-(1)`, `Ar(2)`, `Arp(3)`. The `chem.inp` `ELEMENTS` block is
`Ar` and `E` only; grep for `He`/`Ne`/`N2` across `chem.inp` returns nothing. The prior run saw
all three as non-reactive bath gases (core 6 species). **`1458a4592` works: the argon-only claim
now needs no bath-gas caveat.**

### 2. Does the Cantera `ck2yaml` failure still occur? — **NO (was YES)**

`could not convert string to float: 'e-(1)'` — **0 occurrences** in `stdout.log`, `stderr.log`,
`RMG.log`. Instead the log carries once:

> `Skipping the Chemkin-to-Cantera translation (cantera_from_ck/). This mechanism contains rate
> laws evaluated at the electron temperature … Cantera's ck2yaml implements no TDEP handler …`

`cantera_from_ck/` is created but **empty** — the translation pass never runs. The prior run had it
run and fail (caught, non-fatal). **`ebb95ad96` works.** `cantera2/chem.yaml` (the native writer,
plasma phase intact) is still produced.

### 3. Does the run still die in `check_model()` at `generate_reverse_rate_coefficient` with `ReactionError: Unexpected kinetics type ElectronCollisionPlasma`? — **YES, identically (the important one)**

Manager's reading **verified** at source and by run. `rmgpy/reaction.py:2003` guards the reverse
branch with `if len(self.products) >= 2:`, **not** `self.reversible`; the argon ionisation
`Ar + e- => Arp + e- + e-` has 3 products, so line 2016 calls `generate_reverse_rate_coefficient()`
on a reaction declared `reversible = False`. That code is byte-identical between the two runtimes,
and the wall survived:

```
File "rmgpy/rmg/main.py", line 1390, in execute
    self.check_model()
File "rmgpy/rmg/main.py", line 1916, in check_model
    violator_list = rxn.check_collision_limit_violation(...)
File "rmgpy/reaction.py", line 2016, in ...check_collision_limit_violation
    kr_list.append(self.generate_reverse_rate_coefficient().get_rate_coefficient(...))
File "rmgpy/reaction.py", line 1614, in ...generate_reverse_rate_coefficient
    raise ReactionError("Unexpected kinetics type ...")
rmgpy.exceptions.ReactionError: Unexpected kinetics type
<class 'rmgpy.kinetics.arrhenius.ElectronCollisionPlasma'>; should be one of (...)
```

Same fatal wall, same file/line, same message as the prior run. **This is the answer, not a bug to
fix here** — the fix lands in `rmgpy/reaction.py` (blast radius beyond plasma) and is a separate
ticket.

### 4. Does `MODEL GENERATION COMPLETED` appear? — **NO (same as prior)**

**0 occurrences** in `stdout.log` and `RMG.log`. Model *generation* finished — the deck integrated
to its termination time, core/edge were written to Chemkin and native Cantera — but `check_model()`
(main.py:1390) raises before main.py logs the completion line (main.py:1418, only after the check).
**True python exit code: `PYTHON_EXIT=1`** (captured with `EXIT=$?` immediately after the
interpreter; the background wrapper's trailing `echo` reported "exited with code 0", the same
masking the prior run flagged). Matches the prior run's exit 1.

### 5. Final counts and `chemkin/chem.inp` verbatim

Core **3 species / 1 reaction** (prior: 6 / 1); edge **2 species / 2 reactions** (prior: 2 / 2).
The only count that changed is core species: 6 → 3, exactly the three bath gases removed.

```
ELEMENTS
	Ar
	E
END

SPECIES
    e-(1)
    Ar(2)
    Arp(3)
END

THERM ALL
   300.000  1000.000  5000.000

e-(1)                   E   1               G   298.000  3000.000 1000.00      1
 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00                   4

Ar(2)                   Ar  1               G   200.000  6000.000 1000.00      1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-7.45375000E+02 4.37967000E+00 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-7.45375000E+02 4.37967000E+00                   4

Arp(3)                  Ar  1E  -1          G   100.000  5000.000 1970.48      1
 2.50000241E+00-2.86966665E-14 1.51088780E-17-3.34839323E-21 2.65220781E-25    2
 1.82138651E+05 5.77103097E+00 2.50000241E+00-1.47796363E-15 3.19065367E-18    3
-2.29397207E-21 5.09250534E-25 1.82138651E+05 5.77103097E+00                   4

END

REACTIONS    KCAL/MOLE   MOLES

Ar(2)+e-(1)=>Arp(3)+e-(1)+e-(1)                     2.981e+13 0.597     361.244
    TDEP/e-(1)/   ! ElectronCollisionPlasma exported as a modified-Arrhenius fit of k(Te) over 1.16e+04-1.16e+05 K; the original functional form is not representable in Chemkin

END
```

The prior run's `SPECIES` line read `He   Ne   N2   e-(1)   Ar(2)   Arp(3)`; here it is
`e-(1) Ar(2) Arp(3)`. The reaction line and its `TDEP/` annotation are byte-identical.

### 6. Full model inventory with provenance

| where | species | reaction | provenance |
|---|---|---|---|
| core | `e-(1)`, `Ar(2)`, `Arp(3)` | `Ar + e- => Arp + e- + e-` | deck-declared species + **PlasmaArgon** |
| edge | `[Li](4)`, `[Lip](5)` | `[Li] + e- => [Lip] + e- + e-`; `[Lip] + e- => [Li]` | **PlasmaElectronImpactIonization**, **PlasmaRadiativeRecombination** |

`[Li]`/`[Lip]` **still on the edge** — same as the prior run — from the two lithium libraries the
deck still loads. The only inventory change from the prior run is that the core no longer carries
`He`/`Ne`/`N2`. No species is attributable to RMG-injected bath gas; no air species anywhere
(`PlasmaAir` is not loaded).

## Corrected findings list

The prior report's *Findings* section carried five items. Items 1–4 are the durable run findings;
item 5 is a database provenance note (`PlasmaAir` carries two datasets) — unaffected by runtime and
untouched here (`PlasmaAir` is not even loaded), so it stands unchanged.

| # | prior finding | this run |
|---|---|---|
| 1 | The wall moved to argon itself: `check_model()` → `generate_reverse_rate_coefficient` has no branch for `ElectronCollisionPlasma` | **CONFIRMED.** Identical traceback, same file/line (reaction.py:2016→1614), same exit 1. The reverse-rate/collision-limit check still does not support the plasma rate law. Owner: RMG-Py; separate ticket. |
| 2 | RMG injects N2 (and He, Ne) into every model unconditionally | **STALE-RUNTIME ARTIFACT.** `1458a4592` stopped this for plasma jobs; core is now argon-only. The prior observation was a property of the 28-behind runtime, not of `plasma`. |
| 3 | Keeping the two lithium libraries puts `[Li]`/`[Li+]` on the edge | **CONFIRMED.** `[Li](4)`/`[Lip](5)` on the edge from `PlasmaElectronImpactIonization` + `PlasmaRadiativeRecombination`, inert, not air chemistry — exactly as before. |
| 4 | Cantera export cannot represent the exported plasma rate (`ck2yaml` fails on `TDEP/`) | **STALE-RUNTIME ARTIFACT** (as a *run event*). `ebb95ad96` now skips the Chemkin→Cantera pass for plasma mechanisms, so the failure no longer occurs. The underlying *limitation* it pointed at is real and unchanged — `ck2yaml` still has no `TDEP` handler — but on this runtime RMG no longer invokes it, so nothing fails. |

Net: **two of the four durable findings were stale-runtime artifacts** (findings 2 and 4), each
already fixed on `plasma`; **two survive** (findings 1 and 3). The single remaining hard wall is
finding 1 — the argon reverse-rate check — which is the answer this re-run was built to isolate,
now that the bath-gas and Cantera noise is gone.

## Verification

| Verifier | result |
|---|---|
| 1. `git log --oneline HEAD..plasma \| wc -l` = 0 | **PASS** — output `0`, quoted above |
| 2. Built from source, exit quoted | **PASS** — `REAL_BUILD_EXIT=0`; fresh .so, distinct inode |
| 3. Six questions answered + compared | **PASS** — above |
| 4. Corrected findings list | **PASS** — above |
| 5. `pytest test/rmgpy/reactionTest.py` green, count quoted | **PASS** — `81 passed, 2 skipped in 8.43s`, exit 0 |

## Machine notes

`/proc/loadavg` at build start: `2.82`; at run start: `11.91`. Build ~1 min, run ~2 s wall.
Both build and run captured stdout/stderr to `.log` files (in the scratchpad run dir, not
committed).
