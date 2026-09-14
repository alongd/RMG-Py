# Census of every argon atom type — declared vs measured

Generated from `evidence/probe_census.py` against the loaded
`rmgpy/molecule/atomtype.cpython-39-x86_64-linux-gnu.so` at `f14de2663`, database
`/home/alon/Code/RMG-database-plasma/input`. Raw output: `evidence/census.stdout.log`.

Nothing in this document is read off the declarations. Every "measured" column is the result of
applying the primitive to a concrete `Atom` and asking `get_atomtype` what came out.

---

## 1. The five argon leaves, as declared

All of `all_double`, `r_double`, `o_double`, `s_double`, `triple`, `quadruple`, `benzene` are `[0]`
for every leaf, so only the three columns perception can actually use are shown.

| type | `single` | `lone_pairs` | `charge` | meaning |
|---|---|---|---|---|
| `Ar`   | `[]` (wildcard) | `[]` | `[0,1,2]` | the generic parent; `specific=['Ar0','Ar0s','Ar0e','Ar+','Ar++']` |
| `Ar0`  | `[0]`     | `[4]` | `[0]` | bond-free neutral ground state |
| `Ar0s` | `[1]`     | `[3]` | `[0]` | singly-bonded neutral — the neutral half of Ar₂⁺ |
| `Ar0e` | `[0]`     | `[3]` | `[0]` | **new in I-222** — bond-free neutral, metastable Ar\* |
| `Ar+`  | `[0,1]`   | `[3]` | `[1]` | cation, bare or bonded |
| `Ar++` | `[0,1,2]` | `[3]` | `[2]` | dication |

Every leaf carries `generic=['R','R!H','R!H!Val7','Rx','Rx!H','Ar']` and `specific=[]`. No two can
match one atom: `Ar0` is separated from `Ar0e` by lone pairs, `Ar0s` from `Ar0e` by `single`, and
`Ar+`/`Ar++` from all three neutrals by charge.

---

## 2. Perception is five times wider than construction — for every argon type

Exhaustive sweep of bond-free argon over `u0..u4 × p0..p4 × c-1..+3` (100 states), asking both what
`get_atomtype` returns and whether the adjacency list will build it:

| type | (u,p,c) triples that perceive as it | of those, constructible |
|---|---|---|
| `Ar0`  | 5 | **1** (`u0 p4 c0`) |
| `Ar0e` | 5 | **1** (`u2 p3 c0`) |
| `Ar+`  | 5 | **1** (`u1 p3 c+1`) |
| `Ar++` | 5 | **1** (`u0 p3 c+2`) |

`get_atomtype` ignores radical electrons entirely, so each type answers for all five `u` values at
its own `(p, c)`; the adjacency list's valency check admits exactly one of them.

**This corrects a framing in the main report.** I described the `u`-blind widening as "the one
behaviour I-222 knowingly widens". The ratio is 5:1 for `Ar0`, `Ar+` and `Ar++` *as well*, and has
been since they were declared. `Ar0e` did add four new perceive-but-not-construct states — `u0`,
`u1`, `u3`, `u4` at `p3 c0`, which previously raised `AtomTypeError` — but it did not introduce the
pattern. It is a fifth instance of a family-wide property, not a new kind of hazard.

---

## 3. The action graph: what the table declares vs what the primitive does

Only six of the fifty argon action slots are declared non-empty. Here are all six, each against the
measured result of applying its primitive to every concrete representative of that type:

| # | declared edge | primitive | measured result | verdict |
|---|---|---|---|---|
| 1 | `Ar+.increment_charge = ['Ar++']`   | `GAIN_CHARGE` on `Ar+` | `Ar++` (bare and bonded) | **true** |
| 2 | `Ar++.decrement_charge = ['Ar+']`   | `LOSE_CHARGE` on `Ar++` | `Ar+` (bare and bonded) | **true** |
| 3 | `Ar0.increment_charge = ['Ar+']`    | `GAIN_CHARGE` on `Ar0` | `u0 p4 c+1` → **no type at all** | **FALSE** |
| 4 | `Ar+.decrement_charge = ['Ar0']`    | `LOSE_CHARGE` on `Ar+` | bare → **`Ar0e`**; bonded → **`Ar0s`** | **FALSE** |
| 5 | `Ar0.decrement_lone_pair = ['Ar+']` | `LOSE_PAIR` on `Ar0` | `u0 p3 c+2` → **`Ar++`** | **FALSE** |
| 6 | `Ar+.increment_lone_pair = ['Ar0']` | `GAIN_PAIR` on `Ar+` | `u1 p4 c-1` → **no type at all** | **FALSE** |

**Four of the six declared argon edges are false, and they are two mutually-closing pairs:**

- **the charge pair** — #3 ⇄ #4. Each is the other's inverse, which is why both pass
  `TestActionGraphClosure`: closure checks that the graph is *symmetric*, never that an edge is
  *true*. This is the pair HIGH 3 names.
- **the lone-pair pair** — #5 ⇄ #6. Same structure, same reason it passes, equally false. Not named
  in HIGH 3; found by this census.

The remaining forty-four slots are empty. Many of them are empty while the primitive *does* produce
a valid type — for example `Ar0e` `FORM_BOND` → `Ar0s`, `Ar0s` `BREAK_BOND` → `Ar0e`, `Ar++`
`GAIN_PAIR` → `Ar0`. Those are **omissions**, not false statements: the group path refuses an action
the molecule path performs. A false edge and an omitted edge fail differently, and only the false
ones are in scope here.

> Note on the raw log: the generic `Ar` rows in §3 of `census.stdout.log` are flagged
> `<<< DISAGREES` for all ten actions. That is an artifact of the probe — generic `Ar` has no
> concrete representative in its representative list, so its measured column is empty and compares
> unequal to a non-empty declaration. Generic `Ar` declaring every action self-preserving is
> correct and is the subject of HIGH 1, not a defect.

---

## 4. What corrects each false pair

Every replacement below is the measured result, not a proposal:

**Charge pair.** `LOSE_CHARGE` on `Ar+` has two right answers depending on bonding, and `Ar0` is
neither:

```
Ar+.decrement_charge   ['Ar0']    ->  ['Ar0e', 'Ar0s']
Ar0e.increment_charge  []         ->  ['Ar+']      # GAIN_CHARGE on Ar0e measures Ar+
Ar0s.increment_charge  []         ->  ['Ar+']      # GAIN_CHARGE on both Ar0s variants measures Ar+
Ar0.increment_charge   ['Ar+']    ->  []           # produces an atom with no type
```

**Lone-pair pair.**

```
Ar0.decrement_lone_pair   ['Ar+']  ->  ['Ar++']    # LOSE_PAIR on Ar0 measures Ar++
Ar++.increment_lone_pair  []       ->  ['Ar0']     # GAIN_PAIR on Ar++ measures Ar0
Ar+.increment_lone_pair   ['Ar0']  ->  []          # produces an atom with no type
```

Both sets close both ways by construction, so `TestActionGraphClosure` stays green with no new
allowlist entry.

Risk already checked: **no kinetics family in `RMG-database-plasma` spells any argon leaf.**
`grep -rln "Ar0\|Ar++" input/kinetics/families/` returns nothing, so dropping `Ar0.increment_charge`
cannot orphan an argon-spelled template there. Generic-`R` and generic-`Ar` templates are unaffected
either way, because they never consult a leaf's lists (that is HIGH 1).

---

## 5. ▶ THE QUESTION ◀

I-222's brief made `Ar0`, `Ar0s`, `Ar+`, `Ar++` non-goals. HIGH 3 lifted that **for `Ar+` only** and
asked me to correct `Ar+.decrement_charge` and add a test that the group path and the perception
path agree — with the instruction to stop and ask if agreement needed `Ar0` or `Ar++` as well.

**It does, and more than HIGH 3 anticipated.** Three things make this a decision rather than an
edit:

1. **The charge pair cannot be fixed from `Ar+` alone.** Removing the false `'Ar0'` orphans
   `Ar0.increment_charge`, and the correct targets `Ar0e`/`Ar0s` need their inverses declared back.
   The minimum closed fix touches `Ar+`, `Ar0`, `Ar0s` and `Ar0e`.

2. **The test HIGH 3 asks for forces the second pair too.** A test that "the group path and the
   perception path agree", scoped to the argon family, fails on the lone-pair pair exactly as it
   fails on the charge pair. Fixing only the charge pair means either leaving that test red or
   scoping it down to the single edge HIGH 3 named — and a test narrowed to the one defect already
   known is the kind that misses its neighbour.

3. **Both pairs are the same defect with the same root cause**: an inverse pair that is symmetric
   and therefore invisible to the closure check, while neither direction states what its primitive
   does.

The options, in ascending scope:

| | what changes | agreement test | residual |
|---|---|---|---|
| **A** | both false pairs — 7 assignments across `Ar0`, `Ar0s`, `Ar0e`, `Ar+`, `Ar++` | real set equality, whole argon family | none in the declared edges |
| **B** | charge pair only — 4 assignments across `Ar0`, `Ar0s`, `Ar0e`, `Ar+` | must exclude the lone-pair pair | pair #5/#6 stays false |
| **C** | `Ar+` lift only: `decrement_charge = ['Ar0','Ar0e']` | weakened to "perceived type is among the declared" | keeps a known-false entry; `Ar0s` case still missing |
| **D** | nothing; document only | none | both pairs stay false, now reachable via `Ar0e` |

My recommendation is **A**, on the strength of point 2: it is the smallest change under which the
test actually asked for can be written honestly. Every one of the seven assignments is forced by a
measurement in §3, none of them widens a feature range, and the database suite (118 tests) is the
check that nothing downstream depended on the false edges.

Not in scope under any option: the forty-four omitted edges, generic `Ar`'s self-preserving
declarations, and anything in RMG-database.
