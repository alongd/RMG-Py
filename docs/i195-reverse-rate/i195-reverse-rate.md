# I-195 — the reverse-rate wall on the irreversible argon ionisation

The 5 torr argon deck builds a clean mechanism (core `e-`, `Ar`, `Ar+`, one reaction) and dies at
final validation: `check_model()` computes a **reverse** collision-limit and reverse rate
coefficient for `Ar + e- => Arp + e- + e-`, which is declared `reversible = False`, and
`ElectronCollisionPlasma` has no reverse-rate branch, so `generate_reverse_rate_coefficient()`
raises `ReactionError`. This ticket **measures first**, then takes the fork honestly. **No RMG
source was changed** — see the verdict.

Runtime: built from source at `plasma` tip `3b479a638` (`python setup.py build_ext --inplace -j 8`,
`BUILD_EXIT=0`, fresh `reaction.cpython-39...so`). Database: `RMG-database-i193-plasmaargon` (carries
`PlasmaArgon`). All measurements below ran against this worktree's own rmgpy (confirmed by the run's
`git HEAD = 3b479a638`, Version 4.0.0), forcing it ahead of the shared `easy-install.pth` pin at the
primary checkout via `PYTHONPATH`.

## The four questions, answered by measurement

### Q1 — is `self.reversible` False *at the point `check_model()` sees it*? — **YES**

Instrumented the real run: a `sitecustomize.py` wrapped `RMG.check_model` (pure Python) to log every
core reaction's flags at the exact point the collision-limit loop iterates them
(`main.py:1912`), before delegating to the original. Captured, verbatim, immediately before the crash:

```
Q1-PROBE core rxn='Ar(2) + e-(1) => Arp(3) + e-(1) + e-(1)' reversible=False kinetics=ElectronCollisionPlasma nreact=2 nprod=3
...
File "rmgpy/reaction.py", line 2016, in ...check_collision_limit_violation
File "rmgpy/reaction.py", line 1614, in ...generate_reverse_rate_coefficient
rmgpy.exceptions.ReactionError: Unexpected kinetics type <class 'rmgpy.kinetics.arrhenius.ElectronCollisionPlasma'>
```

The runtime core reaction is `reversible=False`, 2 reactants, 3 products — a library flag confirmed
on the runtime object, not merely in the file. The reverse branch (`reaction.py:2010`, guarded by
`len(self.products) >= 2`, **not** `self.reversible`) is entered because it has 3 products.

### Q2 — is `generate_reverse_rate_coefficient()` meaningless for an irreversible reaction? — **NO. It is a thermodynamic quantity, `k_f/K_eq`, independent of the `reversible` flag.**

`generate_reverse_rate_coefficient()` computes `k_r(T) = k_f(T) / K_eq(T)`, where `K_eq` comes from
species thermochemistry (`reaction.py:1528,1557,...`). It never reads `self.reversible`. Proven on
`reactionTest`'s `reaction2` (acetyl + O2 => acetylperoxy, Arrhenius, full thermo), flipping the flag:

```
reversible=True   kr(1000K)=8.318728e+07   kf/Keq=8.416289e+07   type=Arrhenius
reversible=False  kr(1000K)=8.318728e+07   kf/Keq=8.416289e+07   type=Arrhenius
```

Byte-identical. The `reversible` flag is a **modelling directive** (does the solver integrate the
reverse reaction?), not a physical claim that `k_r = 0`. Microscopic reversibility guarantees
`k_r = k_f/K_eq > 0` at any finite T. A **reverse** collision-limit violation therefore flags that
the specified **forward** rate is thermodynamically inconsistent — it implies an unphysically fast
reverse — which is actionable data-quality information **regardless of the `reversible` flag**.
The reverse check is not meaningless for irreversible reactions.

### Q3 — how many irreversible reactions reach the reverse branch? — **1964** (blast radius)

Measured over **69 kinetics libraries / 20 599 reactions** in `RMG-database-i193-plasmaargon`
(every `libraries/*/reactions.py`, loaded through the real `KineticsDatabase.load_libraries` path):

| bucket | count |
|---|---|
| reactions with `len(products) >= 2` (enter reverse branch) | 17 586 |
| of those, `reversible = False` (**blast radius**) | **1964** |
| ├─ reverse-**supported** kinetics — silently reverse-checked today, no crash | 1941 |
| └─ reverse-**unsupported** kinetics — would crash the run | 23 |
| &nbsp;&nbsp;&nbsp;&nbsp;`TwoTemperaturePlasma` (PlasmaAir) | 19 |
| &nbsp;&nbsp;&nbsp;&nbsp;`ElectronCollisionPlasma` (PlasmaAir 3 + PlasmaArgon 1) | 4 |

So **1941 ordinary irreversible reactions already receive a silent reverse collision-limit check
today** (exactly as the contract predicted for irreversible Arrhenius with ≥2 products), and Q2
proves that check is legitimate. Only **23** plasma-rate-law reactions crash, all confined to
kinetics types with no reverse implementation. (Top contributors to the 1964: CH3Cl 1156,
CurranPentane 369, CF2BrCl 134, Narayanaswamy 93, Chernov 45.)

### Q4 — is the forward check affected by the change under consideration? — **NO**

The forward branch (`reaction.py:2003-2009`, `2017-2022`) uses `calculate_coll_limit(reverse=False)`
and `get_rate_coefficient` — it never calls `generate_reverse_rate_coefficient`. Only the reverse
branch (`2010-2016`, crash at `2016`) does. Any fix scoped to the reverse call leaves the forward
check byte-for-byte unchanged.

## The verdict — (A) does NOT clear its bar

**(A) Guard the reverse branch with `self.reversible`.** Its bar: *show the suppressed reverse
checks on irreversible reactions are meaningless.* **They are not** (Q2): the reverse rate is a
flag-independent thermodynamic quantity, and a reverse violation flags a thermodynamically
inconsistent forward rate. (A) would suppress **1941 legitimate, currently-running** reverse checks
across ordinary combustion chemistry to dodge **23** plasma crashes. On this campaign the burden of
proof is on weakening a refusal; (A) cannot carry it. **(A) fails its bar, and is not implemented.**

**(B) Give `ElectronCollisionPlasma` reverse-rate support** — **GATED to the owner, forbidden, and
not attempted.** It is also the physically hard fix: the forward rate is evaluated at the electron
temperature `Te` while `K_eq` is a gas-temperature quantity, so the reverse of a two-temperature
plasma reaction is *not* simply `k_f/K_eq(T_gas)`. Detailed balance across two temperatures is a
research question — which is exactly why it is the owner's.

**(C) Degrade gracefully** — the collision-limit check should still run, but when a kinetics type has
no reverse implementation it should **warn and skip that one reaction's reverse direction** rather
than raise and abort the whole `check_model()`. This is strictly better than (A): it preserves all
1941 legitimate reverse checks and only skips the 23 genuinely un-computable ones (whose reverse rate
cannot be formed without the gated physics of B). It converts an uncaught crash — which today
produces *no* verdict for *any* reaction once it fires mid-loop — into a completed check plus a
logged warning. **Recommended, but not implemented — see below.**

## Why nothing was implemented

(C) changes *when* the collision-limit check fires (it skips un-computable kinetics), which lands on
the ticket's gate: *"Any change that makes a guard, check or warning fire less often, beyond the
single `self.reversible` condition — Stop and report."* (C) is therefore a **new decision**, not this
ticket, and is surfaced to the manager rather than implemented. (A) fails its bar; (B) is gated.
No defensible fix is inside this ticket's lane, so the honest deliverable is the measurement, the
verdict, and the fork.

## The argon re-run (current, unfixed behaviour)

`python rmg.py input.py`, true exit captured from the interpreter itself:

```
PYTHON_EXIT=1
```

- `MODEL GENERATION COMPLETED`: **0 occurrences** (`stdout.log`, `RMG.log`) — `check_model()` raises
  before the completion line.
- Core: **3 species / 1 reaction** (`e-(1)`, `Ar(2)`, `Arp(3)`; `Ar + e- => Arp + e- + e-`).
- Edge: **2 species / 2 reactions** (`[Li](4)`, `[Lip](5)` — inert lithium, not ours).
- `chemkin/chem.inp`: `SPECIES {e-(1), Ar(2), Arp(3)}`; the one reaction with its `TDEP/e-(1)/`
  annotation. Identical to the I-194 record.

A further wall is not reached — this *is* the wall, unmoved, because no fix was applied.

## Verification

| # | check | result |
|---|---|---|
| 1 | Q1–Q4 answered with measurements | **PASS** — above, each with file:line or command+output |
| 2 | verdict stated with blast-radius number | **PASS** — (A) fails; blast radius 1964 (1941 legit / 23 crash) |
| 3 | (A) implemented + RED/GREEN regression test | **N/A** — conditional on (A) clearing its bar; it does not |
| 4 | argon re-run reported exactly | **PASS** — `PYTHON_EXIT=1`, no completion line, core 3/1, edge 2/2, chem.inp above |
| 5 | `pytest test/rmgpy/reactionTest.py` green, baseline + after | **PASS** — `81 passed, 2 skipped`; no source changed, so after == baseline |
| 6 | one suite per process | **PASS** — reactionTest run alone |

## Durable findings

- `check_collision_limit_violation` checks **both directions** by design (`main.py:1910` comment),
  gated on reactant/product count, **never** on `reversible`. The reverse rate it forms is
  `k_f/K_eq`, a thermodynamic quantity independent of the `reversible` flag (Q2). A reverse
  violation on an *irreversible* reaction is a real signal about the *forward* rate's
  thermodynamic consistency — not noise.
- 1964 irreversible ≥2-product reactions across the database enter the reverse branch; 1941 are
  silently, legitimately reverse-checked today; only 23 (all `TwoTemperaturePlasma` /
  `ElectronCollisionPlasma`, in PlasmaAir + PlasmaArgon) crash for lack of a reverse implementation.
- The physically-correct reverse for a two-temperature plasma rate is **not** `k_f/K_eq(T_gas)`
  (forward at `Te`, `K_eq` at `T_gas`); this is the owner's gated (B) work.
- Shared-env gotcha: `easy-install.pth` pins `import rmgpy` at the **primary** checkout
  `/home/alon/Code/RMG-Py`. A worktree run must lead `PYTHONPATH` with its own path, or an early
  `sitecustomize` import will cache the primary's rmgpy before the script dir is on `sys.path`.
