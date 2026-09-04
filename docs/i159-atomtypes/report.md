# i159 — Repair the seven restored atom types

Branch `i159-atomtypes`. Base commit under repair: `dbd131221` ("Restore the seven atom types
plasma was missing relative to branch 99"), which added `Li-`, `Na-`, `K-`, `Ar0`, `Ar+`, `Ar++`,
`N3dc` to `rmgpy/molecule/atomtype.py`.

## Verdicts

| Type(s) | Verdict | One-line reason |
|---|---|---|
| `Ar0` / `Ar+` / `Ar++` | **KEEP** | Faithful port of branch 99; type the only two argon signatures the database carries; safe (see argon section) |
| `Li-` / `Na-` / `K-` | **DROP** | Removal commit b52045138's central reason — nothing uses an alkali anion — is confirmed still true (zero uses in the database) |
| `N3dc` | **DROP** | Zero database use, absent from mainline, unwired even in source branch 99, off-theme nitrogen, an isolated `ActionError` trap |

Net change vs the base: the file is reduced to the argon three. `git diff 4f18dc389` (the
pre-restore parent) on `atomtype.py` shows only argon additions and `nonSpecifics` losing `Ar`.

---

## 1. N3dc — wiring, or why removed

**Removed.** `N3dc` (a doubly-bonded, charge-separated nitrogen: `single=[0]`, `all_double=[1]`,
`lone_pairs=[1]`, `charge=[+1]`) was registered in `N`'s `specific` list and reachable through
`get_atomtype`, but **all ten of its action lists were empty** — it had no `set_actions` call at
all. Confirmed by probe: every one of `increment_bond … decrement_charge` returned `[]`. Any recipe
that types an atom as `N3dc` and applies any action raises `ActionError`.

Why removed rather than wired:

- **No database consumer.** `grep -rn "N3dc" /home/alon/Code/RMG-database-plasma/input` → **0 hits**.
  Nothing in the carried chemistry produces or names an `N3dc` atom.
- **Absent from mainline.** `git show main:…/atomtype.py | grep -c N3dc` → 0; same for
  `official/main` and `official/master`. This repo lineage has never carried `N3dc`.
- **Unwired even in its source.** `dbd131221` ports `N3dc` from branch 99. Branch 99 *also* has no
  `set_actions` for `N3dc` (`git show 99:…/atomtype.py` shows the declaration, no action line). So
  there is no authoritative in-lineage wiring to copy, and it is dead in 99 too.
- **Off-theme.** Every other member of this change set is an alkali or noble-gas plasma ion. `N3dc`
  is nitrogen; it does not belong to the argon/alkali plasma work.
- **A closure check cannot see this defect.** An all-empty node has no edges, so it is trivially
  "closed" — the both-ways check (§2) reports zero violations for `N3dc`. Its defect is *isolation
  / reachability*, a different class. Wiring it to satisfy a closure check would be building to
  pass the check; dropping the unused, off-theme type is the correct resolution.

Dropping it returns argon-adjacent nitrogen typing to the known-good base state (the current
database loads and its consistency suite passes with this change — see the DB section). No atom in
the database typed as `N3dc` (zero uses), so no species becomes untypeable.

## 2. Alkali anions — does any carried artifact use one?

**No — zero uses, every category.** Exhaustive search of `/home/alon/Code/RMG-database-plasma/input`:

| Search | Pattern | Hits |
|---|---|---|
| Li with negative charge | `Li … c-1` | 0 |
| Na with negative charge | `Na … c-1` | 0 |
| K with negative charge | `K … c-1` | 0 |
| Group atom-type tokens | `Li-` / `Na-` / `K-` as adjacency atom type or `{Li-,…}` | 0 |
| SMILES | `[Li-]` / `[Na-]` / `[K-]` | 0 |

Every alkali atom in the database carries `c0` or `c+1`, never a negative charge. (The many
`Li-O`/`Li-N`/`Li-like` strings are bond-type descriptors, group labels, and astrophysics
isoelectronic-sequence terms in comments — not anion atom types. Every `c-1` token in the tree
belongs to C, N, O, `e`, or `R!H`.)

Removal commit **b52045138** gave three reasons for dropping `Na-`/`K-`: (a) no carried reaction,
group or dictionary uses an alkali anion; (b) this line's lithium had no `Li-` either; (c) as
transcribed they were wrong — the sample builder produced `Na u1 p0 c-1` bonded to an H. Reason
(c) has been **fixed** in `dbd131221` (see §3), but reason (a) — the load-bearing one — **still
holds**: zero uses. Per the ticket, "zero uses means this branch is re-adding three types nothing
needs," and "reverting a reasoned removal needs positive evidence, and its absence is itself the
answer." There is no positive evidence. Dropped.

## 3. What `make_sample_atom` / `make_sample_molecule` produce for each anion

Measured on `dbd131221`'s wiring, before removal:

| Type | `make_sample_atom` | `make_sample_molecule` |
|---|---|---|
| `Li-` | `Li` radical=0 lone_pairs=1 charge=-1 | `1 Li u0 p1 c-1` (bare atom) |
| `Na-` | `Na` radical=0 lone_pairs=1 charge=-1 | `1 Na u0 p1 c-1` (bare atom) |
| `K-`  | `K`  radical=0 lone_pairs=1 charge=-1 | `1 K u0 p1 c-1` (bare atom) |

This is **different** from the nonsense the removal commit described (`Na u1 p0 c-1` bonded to an
H). The current declarations use `lone_pairs=[1]` and `single=[0]`, so the sample builder emits a
clean closed-shell anion (`u0 p1 c-1`) with no spurious bond. **The transcription defect is fixed.**
But a correctly-transcribed atom type that nothing uses is still a type nothing uses — the fix
removes reason (c) for keeping them out, not reason (a). They are dropped on (a).

## 4. Argon `get_atomtype` — before/after, for every database argon signature

Only **two** argon signatures occur in the database (37 neutral + 4 cation, all concrete species;
**no argon group** uses a generic `Ar` atom-type token anywhere):

| Adjacency list | Occurrences | BEFORE (`Ar` in `nonSpecifics`) | AFTER (`Ar` removed) |
|---|---|---|---|
| `Ar u0 p4 c0` (neutral) | 37 | `Ar` (generic) | **`Ar0`** |
| `Ar u1 p3 c+1` (Ar cation) | 4 | `Ar` (generic) | **`Ar+`** |

Measured by toggling `nonSpecifics` at runtime (it is read as a module global inside
`get_atomtype`). Before, both signatures short-circuit to the lumped generic `Ar`; after, each
resolves to its specific type — exactly the neutral/ionised distinction an argon discharge needs,
and what branch 99 asserts. Both concrete database species type correctly under the new scheme.

**The feature table is wider than the evidence (noted, faithful to branch 99).** The kept types
carry bond-count features that admit *bonded* argon cations — `Ar+` has `single=[0,1]` and `Ar++`
has `single=[0,1,2]` (`Ar0` is `single=[0]`, monatomic-only) — whereas the only argon signatures the
database actually carries are monatomic (`Ar u0 p4 c0`, `Ar u1 p3 c+1`). This is an exact port of
branch 99, not a new claim; it is recorded here so a reader does not have to discover that the
feature table permits argon–X bonds no carried species exercises.

## 5. Do argon groups still match? (constructed case)

**Yes.** The database contains no argon groups, so no real group match is at risk — but the ticket
asks for the constructed case, and it holds. A group written with the generic `Ar` atom type
matches a concrete neutral-argon molecule whose atom types as `Ar0`:

    group  "1 Ar u0"           vs  molecule "1 Ar u0 p4 c0"  ->  is_subgraph_isomorphic = True

The generic/specific hierarchy carries it: `Ar0`'s `generic` list includes `'Ar'`, so a group
demanding atom type `Ar` is satisfied by an atom typed `Ar0`. Removing `Ar` from `nonSpecifics`
therefore does not orphan any (hypothetical) generic-`Ar` group. This is a guarded regression test:
`test_generic_argon_group_matches_specific_argon_atom` in `atomtypeSevenTest.py`.

### Is argon being singled out correctly, or do He/Ne need the same?

**Singled out correctly.** Argon is this campaign's immediate target (the 5 torr argon-only
simulation `Ar + e- => Ar+ + 2 e-`). `He` and `Ne` remain lumped (empty `specific`,
`charge=[0,1,2]`, still in `nonSpecifics`) — safe by construction, and giving them specific types
would be a port of types outside this ticket's seven-type scope (a non-goal). He+/Ne+ species do
exist in the carried `PlasmaAir` library, so He/Ne specific typing is a **deliberate deferral, not a
defect** — noted here so the next person does not read the asymmetry as an oversight.

## 6. Verdict on the generic `Ar` self-loops

**Intended port; inert only by database contents, not by construction.** `dbd131221` gives the
generic `Ar` all ten actions as self-loops (`increment_bond=['Ar'], …, decrement_charge=['Ar']`),
where plasma previously had none. This is an **exact, faithful port of branch 99** (`git show
99:…/atomtype.py` line 855 is identical). It is the same pattern every generic parent type carries
(`metal`, `alkali`, `C`, …): a recipe acting on a generically-`Ar`-typed atom stays generically `Ar`.

What it permits that was previously refused: **any recipe action on an atom explicitly typed generic
`Ar` now succeeds where it previously raised `ActionError`.** This path is **not** closed by
`get_atomtype` never returning generic `Ar`. `GroupAtom`'s action methods (`rmgpy/molecule/group.py`
~L362–463) read the atom type's action lists **directly**, so any group or template atom *written*
with the `Ar` atom type absorbs the action regardless of what `get_atomtype` would infer for a
molecule. It is inert today only because **no shipped argon group in RMG-database (plasma) uses a
generic `Ar` atom type** — indeed the database has no argon groups at all (§4/§5). That is a
database-contents premise, not a structural guarantee: a future group written with a bare `Ar` atom
type, or a recipe that types an intermediate as generic `Ar`, would exercise the self-loops. The
self-loops are internally consistent for closure (each action's inverse maps `Ar`->`Ar` back), so
they introduce no closure violation regardless.

---

## Argon safety — the load-bearing question (not escalating)

The `nonSpecifics` removal is **safe**; I am not escalating. Evidence:

- Both concrete database argon signatures type correctly to `Ar0`/`Ar+` (§4).
- The database contains **no argon groups**, so the feature-matching path the removal commit worried
  about is never exercised for argon in the carried chemistry.
- The one constructed group case (generic `Ar` vs `Ar0` molecule) matches via the hierarchy (§5).
- Argon thermo/transport library entries are looked up by molecular identity (graph isomorphism),
  not atom-type group matching, so specific typing (`Ar0`) vs generic (`Ar`) does not change library
  resolution.
- The argon change is a faithful port of branch 99, which asserts `Ar0` for neutral argon.

---

## Closure findings (§2) — pre-existing, left as-is

The both-ways action-graph closure check (`TestActionGraphClosure` in `atomtypeTest.py`) is
generated from `ATOMTYPES` and runs over the whole table. On `dbd131221` it found **43** one-way
edges. Six were the alkali-anion sinks introduced by the restore (each anion declared
`decrement_lone_pair`/`increment_charge` toward its neutral, but `Li0`/`Na0`/`K0` declared neither
inverse back — branch 99 wires those inverses, the restore dropped them). Those six vanish with the
anions. The remaining **37 are pre-existing** (base commit `4f18dc389`) and out of scope per the
ticket's non-goal; they are recorded verbatim in `KNOWN_PREEXISTING` so the test passes today yet
fails the instant a *new* asymmetry appears. They span:

- `Mg0d->Mg+`, `Ca0d->Ca+` (increment_charge with no decrement_charge back)
- carbon: `Cs`/`Csc`/`C2sc`/`Cd`/`CO`/`CS`/`Ct`/`Cbf`/`C2tc` lone-pair and bond edges
- nitrogen: `N0sc->N1sc`, `Nm1`/`Nm2`/`Nm3` self lone-pair increment with no decrement
- oxygen: `Oa`/`O0sc`/`O2d` lone-pair edges
- phosphorus: `P1sc`/`P3t`/`P5s`/`P5d` bond/lone-pair edges
- sulfur: a dozen `S4*`/`S6*` bond, form_bond, lone-pair and charge edges

None involve argon or the alkali metals; the argon and alkali families close both ways after the
fix (asserted by `test_argon_and_alkali_families_close_both_ways`).

---

## Verification

Environment: `conda activate rmg_env`; `python setup.py build_ext --inplace -j 8` after every edit;
`PYTHONPATH` set to this worktree; module path confirmed as
`/home/alon/Code/RMG-Py-i159-atomtypes/rmgpy` before each result.

### Unit tests (`-o addopts=""` to drop the repo's default coverage flags)

`pytest test/rmgpy/molecule/atomtypeTest.py atomtypeSevenTest.py moleculeTest.py groupTest.py`
-> **312 passed, 3 skipped** (baseline 305 passed / 3 skipped).

Delta **+7**, fully accounted for:
- `atomtypeSevenTest.py` reworked from 7 no-crash membership tests to 11 (3 argon sample-atom, 4
  dropped-type absence, 2 concrete-argon typing, 1 group-match, 1 argon closure): **+4**.
- `atomtypeTest.py` gains `TestActionGraphClosure` (3 tests): **+3**.

### Verifier #5 — the reworked `atomtypeSevenTest.py` fails against the tip

Rebuilt with `dbd131221`'s `atomtype.py`, ran the reworked file:
`4 failed, 7 passed` — the four `test_dropped_types_are_absent[Li-|Na-|K-|N3dc]` fail because those
types are still registered on the tip. Restored and rebuilt the fixed `atomtype.py` afterward.

### Database consistency suite — passes on the merge target; 183 figure still unreconciled

`test/database/databaseTest.py` collects **6 tests** (no parametrization, no subtests), so the
ticket's "183 passed" figure does not correspond to this file at pytest granularity. I could not
reconcile the 183 count and report that openly.

On the merge target the suite **passes**. Merging this branch (`af12b4a60`) onto the plasma tip
(`7f052ddf8`), rebuilding, and running `test/database/databaseTest.py` against the current database
(`RMG-database-plasma@plasma c833590ec`): **6 passed, 0 failed in 303.92 s**. The current database
loads fully and every consistency check passes with this atom-type change applied.

**Correction to an earlier draft of this report.** An earlier draft claimed the current database
"cannot load at all" because `PlasmaArgon`'s reaction `Ar + e- => Arp + e- + e-` fails
`Reaction.is_balanced()` (free-electron element-`e` count 1 vs 2), producing 6 errors. **That was
measured against this branch's base engine and is false for the merge target.** This branch's
merge-base with plasma is `4f18dc389`, which predates plasma commit `33350c38d` ("Judge the free
electron by charge, not by element count"). On plasma, `is_balanced()` exempts the free electron
from the per-element census, so all three plasma reactions balance:

    PlasmaElectronImpactIonization | [Li] => [Lip]            | electrons=1  | is_balanced=True
    PlasmaRadiativeRecombination   | [Lip] => [Li]            | electrons=-1 | is_balanced=True
    PlasmaArgon                    | Ar + e- => Arp + e- + e- | electrons=0  | is_balanced=True

The apparent load failure was an artifact of the branch base lagging plasma — **not** a database
defect, **not** a load blocker on the merge target, and nothing to escalate. The `053d48080`
pre-`PlasmaArgon` substitute run that the earlier draft relied on for its "database-neutral" proof
is therefore moot and withdrawn; the direct 6/0 on the current DB replaces it.

The change is **database-neutral**: the consistency suite passes 6/0 with the current engine, the
current database, and this atom-type change applied. Atom typing does not enter the balance
question at all — element symbol and charge come from the adjacency-list tokens, not the atom type.

## What I could not determine

- **The "183 passed" baseline.** `databaseTest.py` is 6 tests; I could not reconcile 183 with any
  collection of it in this worktree. The substantive question it protects — does this change regress
  database consistency — is answered (no): the suite passes **6/0** on the merge target with this
  change applied. Only the exact 183 figure is not reproduced.
