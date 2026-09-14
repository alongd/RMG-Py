# I-222 — An atom type for bond-free neutral argon (metastable Ar\*)

Branch `i222-metastable-argon-atomtype`, worktree
`/home/alon/Code/RMG-Py-i222-metastable-argon-atomtype`, base `78f306665`.
Nothing pushed, nothing merged.

**Outcome: the type is added, and the argon action graph is repaired.** `Ar u2 p3 c0` — metastable
argon — now builds as a `Molecule` and as a `Species` through the ordinary path and types as
**`Ar0e`**. The tolerant path no longer degrades it to the generic wildcard `R`. Four of the six
declared argon action edges were measurably false and are now corrected.

> ## Scope notice — this branch changes pre-existing argon behaviour
>
> **It is not only an addition.** Alongside the new `Ar0e`, it rewrites the action table of all four
> argon leaves that already existed. Anything already written against `Ar0`, `Ar0s`, `Ar+` or `Ar++`
> — any group template advanced through the group path by a charge or lone-pair action — is affected
> independently of metastable argon.
>
> | leaf | before | after |
> |---|---|---|
> | `Ar0` | `decrement_lone_pair=['Ar+']`, `increment_charge=['Ar+']` | `decrement_lone_pair=['Ar++']`, `increment_charge=[]` |
> | `Ar0s` | `increment_charge=[]` | `increment_charge=['Ar+']` |
> | `Ar+` | `increment_lone_pair=['Ar0']`, `decrement_charge=['Ar0']` | `increment_lone_pair=[]`, `decrement_charge=['Ar0e','Ar0s']` |
> | `Ar++` | `increment_lone_pair=[]` | `increment_lone_pair=['Ar0']` |
>
> **Every old edge was false**, so these are corrections, not changes of intent. The arithmetic, one
> line each — `Ar0` is `p4 c0`, bond-free, radical-free:
>
> - **`Ar0` `LOSE_PAIR`**: removing a lone pair leaves six valence electrons paired, so charge
>   recomputes to **+2**. That is `Ar++`. The old edge was off by a charge unit.
> - **`Ar0` `GAIN_CHARGE`**: would be `p4 c+1`, which **no argon leaf admits**. Removing the edge is
>   the only truthful option; ionising ground-state argon is a *compound* of `GAIN_CHARGE` with
>   `LOSE_PAIR`/`GAIN_RADICAL`, and no single primitive stands for it.
> - **`Ar+` `GAIN_PAIR`**: likewise `p4 c-1`, owned by no leaf.
> - **`Ar+` `LOSE_CHARGE`**: leaves `p3 c0`, which is `Ar0e` or `Ar0s` — never `Ar0`, which sits at
>   `p4`. This was the finding an earlier adversarial round raised against the ticket.
>
> The non-goal forbidding edits to `Ar0`/`Ar0s`/`Ar+`/`Ar++` was lifted by the owner for exactly this
> repair, after the census was presented (§3.2). Re-measured in
> `evidence/probe_spar48.stdout.log`; pinned by
> `TestArgonActionPathsAgree::test_no_declared_argon_edge_is_false`.

This document covers **three rounds**. The second was opened by an adversarial review that found
three things, all of which held up:

| | finding | where |
|---|---|---|
| HIGH 1 | `Ar0e` is **not** inert. Empty `set_actions` binds the group path only; generic `Ar`/`R` templates match metastable argon and the molecule path acts on it regardless. My first-round claim to the contrary was wrong. | §3.1 |
| HIGH 2 | An alkaline-earth family can now reach metastable argon. **Not repaired here** — a database change, gated, taken by the owner; recorded as a consequence of this change. | §8.1 |
| HIGH 3 | `Ar+.decrement_charge` named a type its primitive cannot produce. The census found a second, equally false pair. Both repaired under a lift of the original non-goal. | §3.2 |

Two claims from the first round are **corrected** rather than quietly edited: the inertness argument
(§3) and "the one behaviour I-222 knowingly widens" (§8). A third, a `§8`/`§9` contradiction about
whether wider suites had been run, is resolved in favour of §9 — they were run, twice.

The third round found **no correctness defect**. It returned four disclosure items — three of them
about what this branch *says* versus what it *does* — and changed no behaviour:

| | item | where |
|---|---|---|
| 1 | The five-leaf action-table rewrite was announced only in a commit body, where a reviewer would have to diff for it. | the scope notice above |
| 2 | `Ar+.decrement_charge = ['Ar0e','Ar0s']` states a bond-conditioned result as an unconditional union, and the group path both overbroadens and unorders it. | §3.4 |
| 3 | `Ar0e` is a metastability label that cannot constrain `u`, and is now a label a database author can type. | §2, and `Ar0e`'s own declaration comment |
| 4 | `ARGON_DB_SIGNATURES` is a hand-written literal whose docstring claimed a database-coverage reach it does not have. | §7.1 |

Supporting document: [`argon-atom-type-census.md`](argon-atom-type-census.md), a complete
declared-vs-measured census of all five argon leaves.

One thing in the brief was wrong and is reported first, because it changed what was authorised.

---

## 0. The contradiction in the brief, and how it was resolved

The brief's Fixed constraints said:

> `TestArgonSingleBondNarrowing` (five tests, around `atomtypeTest.py:1397`) **stays green**. If your
> change moves either, stop and report.

Two of those five tests assert *exactly* the behaviour the same brief's Verifier items 3 and 4
require to change. Measured, not predicted — with `Ar0e` declared and the extension rebuilt:

```
FAILED atomtypeTest.py::TestArgonSingleBondNarrowing::test_bond_free_triplet_argon_has_no_atom_type
    atomtypeTest.py:1432: Failed: DID NOT RAISE <class 'rmgpy.exceptions.AtomTypeError'>
FAILED atomtypeTest.py::TestArgonSingleBondNarrowing::test_untypeable_argon_degrades_to_generic_R_when_typing_is_tolerant
    atomtypeTest.py:1454: AssertionError: assert 'Ar0e' == 'R'
```

Work stopped there and the conflict went to the owner rather than being resolved locally.

**These two tests were not obstacles that were removed. They were a decision point that was
reached, and the behaviour they pinned was changed deliberately.** The second test says so itself,
in its own docstring as written for I-218:

> Records a HAZARD, not a desideratum. `Molecule.update_atomtypes` catches AtomTypeError and assigns
> `ATOMTYPES['R']` when `raise_exception` is False, with the logging gated behind a separate
> `log_species` flag — so an untypeable atom can silently become the wildcard that matches
> everything. `Species.get_resonance_hybrid` passes both flags off. **This test exists so that
> whoever adds an argon metastable type sees the consequence in CI rather than in a generated
> model; if the silent fallback is ever fixed, change it.**

I-222 is that ticket and `Ar0e` is that type. Both tests are tripwires built to be tripped by
exactly this change. The owner verified the account against the file before authorising, and the
authorised scope was precisely: rewrite those two in place to the post-`Ar0e` truth, rename them so
their names no longer assert the old behaviour, record ticket and commit in each docstring, add
`Ar0e` to `EXPECTED_FAILING_ATOMTYPES` for the measured `make_sample_atom` reason, and touch no
other test.

The other three tests in the class pin the narrowing itself rather than the absence of a metastable
type — `Ar0s.single == [1]`, Ar2⁺ types as `Ar0s`/`Ar+`, `Ar0` keeps the bond-free closed shell.
They are untouched and green. The class still holds exactly five tests. `Ar0s` was **not** widened
and is still in `EXPECTED_FAILING_ATOMTYPES`.

The coverage the second test used to provide was not dropped: the generic-`R` fallback is still
pinned, on an argon state that genuinely has no type
(`TestMetastableArgonAtomType::test_untypeable_argon_still_degrades_to_the_wildcard`, the argon
anion). That test passes both before and against the change — see §6.

---

## 1. The measurement, reproduced

Run before any edit, in `rmg_env`, against base `78f306665`
(`docs/i222-metastable-argon-atomtype/evidence/probe_baseline.py`, output in
`evidence/baseline.stdout.log`):

```
database.directory = /home/alon/Code/RMG-database-plasma/input
loaded atomtype module __file__ = /home/alon/Code/RMG-Py-i222-metastable-argon-atomtype/rmgpy/molecule/atomtype.cpython-39-x86_64-linux-gnu.so
loaded molecule module __file__ = /home/alon/Code/RMG-Py-i222-metastable-argon-atomtype/rmgpy/molecule/molecule.cpython-39-x86_64-linux-gnu.so

Ar0s.single       = [1]
Ar0s.lone_pairs   = [3]
Ar0.lone_pairs    = [4]
existing Ar types = ['Ar', 'Ar0', 'Ar0s', 'Ar+', 'Ar++']

Molecule().from_adjacency_list('1 Ar u2 p3 c0') -> RAISES AtomTypeError
    Unable to determine atom type for atom Ar.., which has 0 single bonds, 0 double bonds (0 to O,
    0 to S, 0 others), 0 triple bonds, 0 quadruple bonds, 0 benzene bonds, 3 lone pairs, and +0 charge.
Species().from_adjacency_list('1 Ar u2 p3 c0') -> RAISES AtomTypeError
    (same message)

tolerant update_atomtypes(log_species=False, raise_exception=False) -> atomtype = R
  and that argon matches "1 *1 R u[1,2,3,4] px c[0,+1,...]" -> True
```

**Every line of the brief's block reproduces.** The only difference is cosmetic: the exception text
reads `+0 charge`, not `0 charge`.

The four-suite baseline also reproduces exactly: `atomtypeTest` **46**, `moleculeTest` **195**,
`groupTest` **69**, `atomtypeSevenTest` **11**; together **318 passed / 3 skipped**
(`evidence/baseline_collect.stdout.log`, `evidence/baseline_suite.stdout.log`).

A fresh worktree carried zero `.so` files and no `rmgrc`. Both were created before measuring:
`cp rmgrc.template rmgrc` (pointing at `../RMG-database-plasma/input`, printed as the first line of
every log) and a full `python setup.py build_ext --inplace`.

---

## 2. The declaration

`rmgpy/molecule/atomtype.py`:

```python
ATOMTYPES['Ar0e'] = AtomType('Ar0e', generic=['R', 'R!H', 'R!H!Val7', 'Rx', 'Rx!H', 'Ar'], specific=[],
                            single=[0], all_double=[0], r_double=[0], o_double=[0], s_double=[0],
                            triple=[0], quadruple=[0], benzene=[0], lone_pairs=[3], charge=[0])
```

Zero bonds of every order, three lone pairs, neutral — as specified. `'Ar0e'` was also added to
`ATOMTYPES['Ar'].specific` (without which `get_atomtype` can never return it: it only ever iterates
`ATOMTYPES[atom_symbol].specific`) and to the `specific` lists of `R`, `R!H`, `R!H!Val7`, `Rx` and
`Rx!H`, matching every other argon leaf. Both directions of every link are asserted by
`test_ar0e_is_linked_into_the_type_hierarchy_both_ways`.

### The name

**`Ar0e`** — argon, charge 0, electronically excited.

- The file's stated convention is `<element> <valence> <characteristic bonds> <charge(optional)>`,
  but the argon family does not follow it: the `Ar0s` comment says outright that "0 is the charge,
  s the single bond", following `Mg0s`/`Ca0s`. `Ar0e` is consistent with the family it joins.
- `e` is not a bond descriptor. Every trailing letter RMG uses for bonds is one of `s d t q b`
  (plus `c` for charged, `a` for atomic), so `e` cannot be misread as a bond that `Ar0e`, which has
  none, would then be claiming.
- **`Ar2` was considered and rejected.** Under the documented *valence* reading it is the correct
  name — the atom has exactly two electrons free to remain as radicals — and it would have been my
  choice in another file. But this file, and this campaign, discuss the Ar₂⁺ dimer constantly
  (`Ar0s`–`Ar+` is that dimer), and an atom type called `Ar2` would be read as the dimer by every
  later reader. Names that are right and misleading lose to names that are arguable and clear.
- **`Ar*` was rejected** on the constraint the brief names: the label is a dictionary key written
  verbatim into group adjacency lists, and `*` is the adjacency-list grammar's atom-label sigil
  (`1 *1 R u0 ...`). It parses today only because labels are detected by a *leading* `*`; that is
  too thin a margin for a permanent key. `Ar0e` is alphanumeric.
- **`Ar0p3` was considered** — it names the actual perception discriminator — and rejected as
  unlike anything else in the table. The point it would have made is made in the code comment
  instead.

### Mutual exclusivity

`get_atomtype` returns the **first** entry of `ATOMTYPES['Ar'].specific` whose every feature range
contains the atom's value, matching on bond counts, lone pairs and charge only. So two argon types
are safe iff no atom can satisfy both — which is also what makes the order of that list unable to
change any answer:

| type   | `single`    | `lone_pairs` | `charge` | separated from `Ar0e` by |
|--------|-------------|--------------|----------|--------------------------|
| `Ar0`  | `[0]`       | `[4]`        | `[0]`    | `lone_pairs` (4 vs 3)    |
| `Ar0s` | `[1]`       | `[3]`        | `[0]`    | `single` (1 vs 0)        |
| `Ar0e` | `[0]`       | `[3]`        | `[0]`    | —                        |
| `Ar+`  | `[0,1]`     | `[3]`        | `[1]`    | `charge` (+1 vs 0)       |
| `Ar++` | `[0,1,2]`   | `[3]`        | `[2]`    | `charge` (+2 vs 0)       |

All ten pairs, not just the four against `Ar0e`, are checked by
`test_no_two_argon_types_can_match_one_atom`, which is generated from the declarations rather than
written out, so it also covers pairs a later edit creates. `Ar0s` is what makes the `single` column
load-bearing: it is `Ar0e` with one bond, and the two are separated by nothing else. Widening
`Ar0s` back to `single=[0,1]` — the change I-218 exists to prevent — would now collide with `Ar0e`
directly, and that test fails on it.

The argument is also driven through `get_atomtype` on the concrete atoms
(`test_each_realizable_argon_atom_types_specifically`): `Ar0` `u0 p4 c0`, `Ar0e` `u2 p3 c0`, `Ar+`
`u1 p3 c+1`, `Ar++` `u0 p3 c+2`, and `Ar0s` as the bonded half of Ar₂⁺.

### Which `u` states the type admits — measured, not asserted

`get_atomtype` ignores radical electrons entirely, so `Ar0e` cannot discriminate by `u`. What
refuses the other states is the adjacency list's valency check. Measured through that path
(`evidence/actions.stdout.log`):

```
'1 Ar u0 p3 c0' -> InvalidAdjacencyListError: Invalid valency for atom Ar (Ar0e) ...
'1 Ar u1 p3 c0' -> InvalidAdjacencyListError: Invalid valency for atom Ar (Ar0e) ...
'1 Ar u2 p3 c0' -> OK, atomtype = Ar0e
'1 Ar u3 p3 c0' -> InvalidAdjacencyListError: Invalid valency for atom Ar (Ar0e) ...
'1 Ar u4 p3 c0' -> InvalidAdjacencyListError: Invalid valency for atom Ar (Ar0e) ...
```

Exactly one radical state is constructible, `u2` — the triplet metastable. That is charge balance,
not the declaration: neutral Ar brings 8 valence electrons, three lone pairs consume six, no bond
consumes any, so two are left unpaired.

Asked the other way round — by **direct perception** rather than by adjacency-list rejection —
`get_atomtype` returns `Ar0e` for all five `u` values at `p3 c0`
(`test_perception_admits_unrealizable_u_states`). The two measurements are the same fact from
either side: perception admits five, construction admits one. `Ar0`, `Ar+` and `Ar++` have the
identical 5:1 ratio (§8).

**So the label cannot mean "metastable" by itself, and it is now a label database authors can
type.** That `get_atomtype` ignores `u` was known from the brief; what making `Ar0e` a real
dictionary key adds is that a group adjacency list may now spell it. Group adjacency lists get no
valency check, so `1 Ar0e ux p3 c0` is legal, and matches `u1` argon — ordinary, non-metastable —
as readily as the `u2` it appears to name. Measured at `p3 c0`: that spelling matches `u1` and `u2`;
`u0` renormalises to `p4` and types `Ar0`; `u3` and `u4` fall to generic `R`
(`evidence/probe_spar48.stdout.log` §3, `test_group_spelled_ar0e_matches_more_than_the_metastable_triplet`).

A group meaning metastable argon must therefore write `u2` explicitly — `1 Ar0e u2 p3 c0` — and the
warning saying so lives in `Ar0e`'s own declaration comment in `atomtype.py`, where an author
looking the label up will meet it, not only here.

---

## 3. `set_actions` — and what it does *not* control

> **This section was rewritten in the rework round.** Its first version argued that `Ar0e` should
> declare no action edges, and drew from that the implication that `Ar0e` was therefore inert. The
> first half was under-measured and the second half was simply wrong. Both are corrected below; the
> full argon census is in [`argon-atom-type-census.md`](argon-atom-type-census.md).

### 3.1 Two paths apply an action, and only one reads these lists

- **the GROUP path** — `GroupAtom.apply_action` maps a group atom's types through `set_actions` and
  raises `ActionError` on an empty list. This is what generates and extends templates.
- **the MOLECULE path** — `Atom.apply_action` mutates `radical_electrons` / `lone_pairs` / `charge`
  / bonds directly, and the result is re-typed afterwards by `get_atomtype`. It **never** consults
  `set_actions`.

So an empty list stops a template *spelled with that leaf* from being advanced through the group
graph, and stops nothing at all on a concrete molecule. Measured, both halves
(`test_empty_set_actions_binds_the_group_path_only`):

```
GroupAtom([Ar0e]).apply_action(FORM_BOND)  -> ActionError: Unknown atom type produced from set ...
ReactionRecipe([FORM_BOND]) on two concrete Ar0e atoms -> succeeds; both atoms retype to R
```

**`Ar0e` is therefore not inert, and empty lists never made it so.** Generic `Ar` (`atomtype.py:884`)
and generic `R` declare *every* action self-preserving — `increment_bond=['Ar']`, `form_bond=['Ar']`,
`increment_charge=['Ar']`, and so on for all ten. Every generic group tested matches a concrete
metastable argon:

```
'1 *1 R ux px cx'        -> matches Ar0e: True
'1 *1 R u[2,3,4] px cx'  -> matches Ar0e: True
'1 *1 Ar ux px cx'       -> matches Ar0e: True
'1 *1 Rx ux px cx'       -> matches Ar0e: True
'1 *1 R!H ux px cx'      -> matches Ar0e: True
```

and the recipe then acts on the concrete atom through the molecule path, which the leaf's own lists
cannot veto. `test_a_generic_template_matches_and_reacts_metastable_argon` pins this end to end.

### 3.2 The declared edges must state what the primitive produces

`TestActionGraphClosure` checks that the graph is **symmetric** — every edge has its inverse. It
cannot check that an edge is **true**. Four of the six declared argon edges were false in exactly
the way that hides from a symmetry check: as two mutually-closing pairs, each entry the inverse of
the other, neither matching the primitive. The census measured all six:

| declared edge | primitive measures | verdict |
|---|---|---|
| `Ar+.increment_charge = ['Ar++']` | `Ar++` | true |
| `Ar++.decrement_charge = ['Ar+']` | `Ar+` | true |
| `Ar0.increment_charge = ['Ar+']` | **no type at all** | FALSE |
| `Ar+.decrement_charge = ['Ar0']` | bare → **`Ar0e`**, bonded → **`Ar0s`** | FALSE |
| `Ar0.decrement_lone_pair = ['Ar+']` | **`Ar++`** | FALSE |
| `Ar+.increment_lone_pair = ['Ar0']` | **no type at all** | FALSE |

The owner lifted the `Ar0`/`Ar0s`/`Ar+`/`Ar++` non-goal for this repair after the census was
presented. Seven assignments, each the measured result:

```python
Ar+.decrement_charge      ['Ar0']  ->  ['Ar0e', 'Ar0s']   # two answers; bonding decides which
Ar0e.increment_charge     []       ->  ['Ar+']            # metastable argon ionising
Ar0s.increment_charge     []       ->  ['Ar+']            # the bonded half of Ar2+ ionising
Ar0.increment_charge      ['Ar+']  ->  []                 # produces an atom with no type
Ar0.decrement_lone_pair   ['Ar+']  ->  ['Ar++']           # a bare pair loss costs two charges
Ar++.increment_lone_pair  []       ->  ['Ar0']            # its inverse
Ar+.increment_lone_pair   ['Ar0']  ->  []                 # produces an atom with no type
```

Both sets close both ways by construction, so every closure test stays green with no new allowlist
entry. `TestArgonActionPathsAgree::test_no_declared_argon_edge_is_false` is the assertion closure
could not make: for every argon leaf and every declared edge, the declared target set must equal the
set of types the primitive actually produces from a concrete atom of that type.

`Ar0e` ends with exactly one declared edge, `increment_charge -> Ar+`. `GAIN_PAIR`, `LOSE_PAIR` and
`LOSE_CHARGE` land on states no argon type owns. `FORM_BOND` does give `Ar0s`, and is left
undeclared: its inverse would be a bond-order edge on a noble gas, which is chemistry this branch
has no consumer for. The radical self-edges are omitted as they are for every other argon leaf.

### 3.3 Still not repaired

The **forty-four empty slots** where the primitive does produce a valid type — `Ar0e` `FORM_BOND` →
`Ar0s`, `Ar0s` `BREAK_BOND` → `Ar0e`, `Ar++` `GAIN_PAIR` → `Ar0` and others. Those are *omissions*,
not false statements: the group path refuses an action the molecule path performs. A false edge and
an omitted edge fail differently, and only the false ones were in scope. Generic `Ar`'s
self-preserving declarations are likewise untouched.

### 3.4 `Ar+.decrement_charge` is a union standing for a condition the grammar cannot express

The repaired edge names **two** targets, and which one is right depends on the bond count:

```
bare Ar+            (u1 p3 c+1, no bonds)      --LOSE_CHARGE-->  Ar0e
singly-bonded Ar+   (u0 p3 c+1, one single)    --LOSE_CHARGE-->  Ar0s
```

An action list is a set of labels with no grammar for a condition, so both are named. As a *set*
that is correct, and it is the only truthful thing the table can say.

**On the molecule path this costs nothing.** The recipe mutates the concrete atom and
`get_atomtype` re-perceives it, arriving at whichever answer that atom's bond count earns. Nothing
consults the list (§3.1).

**The residual is on the group path, and it is twofold:**

- **Overbreadth.** `GroupAtom._lose_charge` (`group.py:386`) maps the group atom through this list
  and keeps *both* entries — it never inspects bond count. A group derived by applying
  `LOSE_CHARGE` to an `Ar+` group atom therefore means "`Ar0e` or `Ar0s`" where the molecule it
  stands for meant one of them. Measured: the resulting atom-type list is `['Ar0e', 'Ar0s']`.
- **Order-sensitivity.** That method ends `self.atomtype = list(set(atomtype))` (`group.py:408`),
  and `AtomType` defines no `__hash__`, so the set is keyed on object identity and carries **no
  order guarantee**. Two consumers then read element `[0]` and nothing else:
  `GroupAtom.make_sample_atom` (`group.py:829`) and `Group.pick_wildcards`
  (`group.py:2858`, `2863-2864`).

  > **Measured, and weaker than it sounds.** Six fresh interpreters returned `['Ar0e', 'Ar0s']` six
  > times out of six (`evidence/probe_spar48.stdout.log`, §2b/2c). So this is a **missing guarantee,
  > not an observed flip** — which is also why nothing has caught it. The declared order happens to
  > survive; nothing promises it will, and adding a third target or reordering these two may change
  > what those two consumers build.

**Deliberately not repaired here.** Conditioning an action edge on bond count is a change to the
action grammar, well beyond this ticket. A sibling ticket is at present resolving a bug in another
family caused by exactly this class — an atom-type list whose first element decided a sample — so
the note is not hypothetical. Recorded in the table's own comment in `atomtype.py` as well as here,
because the next person to edit that line will be reading the code, not this report.

---

## 4. What changed, by file

`rmgpy/molecule/atomtype.py`
- new `ATOMTYPES['Ar0e']` declaration, with the reasoning above as comments;
- `'Ar0e'` added to `ATOMTYPES['Ar'].specific` and to the `specific` lists of `R`, `R!H`,
  `R!H!Val7`, `Rx`, `Rx!H`;
- a comment block above the argon `set_actions` stating the group-path/molecule-path split (§3.1)
  and the requirement that an edge state what its primitive produces (§3.2);
- the seven action-edge assignments of §3.2, touching `Ar0`, `Ar0s`, `Ar0e`, `Ar+` and `Ar++`. No
  feature range (`single`, `lone_pairs`, `charge`, bond counts) of any pre-existing type is
  changed — `Ar0s` is still `single=[1]`;
- *(third round, comments only)* the `u`-cannot-be-constrained warning in `Ar0e`'s declaration
  comment (§2) and the bond-conditioned-union note above `Ar+`'s `set_actions` line (§3.4). No
  executable line changed; the extension was rebuilt and re-value-asserted so source and `.so`
  stay in step.

`test/rmgpy/molecule/atomtypeTest.py`
- `EXPECTED_FAILING_ATOMTYPES` gains `"Ar0e"` (see §5). `"Ar0s"` stays.
- the two tripwires rewritten in place and renamed, per §0:
  `test_bond_free_triplet_argon_has_no_atom_type` → `test_bond_free_triplet_argon_types_as_ar0e_not_ar0s`,
  `test_untypeable_argon_degrades_to_generic_R_when_typing_is_tolerant` →
  `test_metastable_argon_no_longer_degrades_to_generic_R_when_typing_is_tolerant`;
- new class `TestMetastableArgonAtomType`, 15 tests;
- new class `TestArgonActionPathsAgree`, 8 tests (5 of them one parametrisation).

`test/rmgpy/molecule/atomtypeSevenTest.py` — brought up to the five-leaf reality: the closure test
covered only `Ar`/`Ar0`/`Ar+`/`Ar++`, so `Ar0s` and `Ar0e` sat outside it exactly while this ticket
rewired the charge edges through them; a registry census now fails if a sixth leaf appears; and the
three leaves with no database species are named outright. The third round narrowed
`ARGON_DB_SIGNATURES`' comment and two docstrings to what they actually check — no test added, no
test removed. Details in §7.1.

Nothing under `rmgpy/data/`, no reactor, no `electron_placement.py`, no RMG-database, no family, no
species, no cross-section. The change did not turn out to require any of them.

Every new assertion names the atom type it expects. None of them asserts merely that nothing raised
— that form passes on generic `R`, which is the failure this ticket removes.

---

## 5. `Ar0e` in `EXPECTED_FAILING_ATOMTYPES`

`TestAtomType::test_make_sample_molecule` failed on `Ar0e` for the identical, already-documented
reason `Ar0s` is exempt, and the reason was measured, not assumed:

```
atomtypeTest.py:174: AssertionError: Couldn't make sample molecules for types Ar0e
```

`GroupAtom.make_sample_atom` (`rmgpy/molecule/group.py:886`) takes the first entry of each feature
list and has **no rule for choosing `u`** — it falls back to `default_atom.radical_electrons`, i.e.
0. So for `Ar0e` it builds argon at `p3 u0`, whose `update_charge` gives `c+2` against the declared
`c0`. Balancing a bond-free neutral argon at `p3` requires `u2`, and the sample builder cannot pick
it. `Ar0s` sits in the list for the same defect, recorded there for I-218.

This **adds** an entry; `Ar0s` was not removed, and the `_wip` tests that iterate the list are
unchanged. Teaching `make_sample_atom` to derive `u` from charge balance would fix both at once and
is a real improvement — it belongs to whoever owns `group.py`, not to this ticket.

---

## 6. Revert-and-rerun: the new tests confirmed RED first

Every arm reverts **only** `rmgpy/molecule/atomtype.py`, leaves the tests in place, rebuilds the
extension, and checks the loaded module by value rather than by mtime.

### Arm A — the original round: base `78f306665`, no `Ar0e` at all

```
LOADED: .../rmgpy/molecule/atomtype.cpython-39-x86_64-linux-gnu.so
Ar0e in ATOMTYPES: False
```

`evidence/revert_red.stdout.log`: **13 failed, 43 passed, 2 skipped** — the 11 `Ar0e` assertions
then in the new class, plus both rewritten tripwires.

### Arm B — the rework: previous commit `f14de2663`, `Ar0e` present, old action graph

```
Ar0e present: True | Ar+.decrement_charge = ['Ar0']
```

`evidence/red_arm_a.stdout.log`: **4 failed, 82 passed, 2 skipped**. Exactly the repair assertions
fail, and nothing else does:

```
FAILED TestMetastableArgonAtomType::test_ar0e_declares_exactly_the_ionisation_edge
FAILED TestArgonActionPathsAgree::test_no_declared_argon_edge_is_false[Ar0]
FAILED TestArgonActionPathsAgree::test_no_declared_argon_edge_is_false[Ar+]
FAILED TestArgonActionPathsAgree::test_the_two_repaired_pairs_state_their_measured_targets
```

This is the arm that matters for the rework: it isolates the seven action-edge assignments from the
declaration, and confirms the agreement test fails on precisely the two leaves that carried false
edges (`Ar0` and `Ar+`) while passing on the three that did not.

### Arm C — the whole rework against base `78f306665`

`evidence/red_arm_b.stdout.log`: **27 failed, 59 passed, 2 skipped** across both test files.

### The tests that pass on both sides, by design

Four parametrisations pass in arm B *and* after the change, and that is their point — they are
`Ar0s`-only cases that do not depend on `Ar0e`
(`test_no_declared_argon_edge_is_false[Ar0s]`, `test_added_argon_types_make_correct_sample_atom[Ar0s]`,
`test_added_argon_types_cannot_make_a_sample_molecule[Ar0s]`,
`test_argon_leaves_without_a_database_species[Ar0s]`). So does
`test_untypeable_argon_still_degrades_to_the_wildcard`, the preserved generic-`R` coverage.

**Stated plainly rather than papered over:** `test_empty_set_actions_binds_the_group_path_only` and
`test_a_generic_template_matches_and_reacts_metastable_argon` are behaviour *pins*, not repairs.
They go red in arm C only because `Ar0e` does not exist there, not because the behaviour they
describe changes. A pin of existing behaviour has no tree against which it can legitimately go red,
and manufacturing one would have proved nothing.

After every arm, `atomtype.py` was restored and rebuilt and the suites returned to green.

---

## 7. Suite counts, per file against its own collected total

| file | collected (base `78f306665`) | collected (first round) | collected (now) | result (now) |
|---|---|---|---|---|
| `atomtypeTest.py`      | 46  | 58  | **68** | 66 passed, 2 skipped |
| `moleculeTest.py`      | 195 | 195 | **195** | 194 passed, 1 skipped |
| `groupTest.py`         | 69  | 69  | **69** | 69 passed |
| `atomtypeSevenTest.py` | 11  | 11  | **20** | 20 passed |
| four together          | 321 (318 passed, 3 skipped) | 333 (330 passed, 3 skipped) | **352** | **349 passed, 3 skipped** |

`atomtypeTest.py` grows by 22 over base — 15 in `TestMetastableArgonAtomType` and 8 in
`TestArgonActionPathsAgree`, less one test that the rework merged away. `atomtypeSevenTest.py` grows
by 9 (§7.1). No other file's collected total moves, and no test anywhere goes from passing to
failing. The two skips in `atomtypeTest.py` are the pre-existing `@pytest.mark.skip(reason="WIP")`
sample-molecule tests; the third skip is
`moleculeTest.py::test_count_internal_rotors_dimethyl_acetylene`, also pre-existing.

Both streams were captured for every measurement; `evidence/` holds the matching `*.stdout.log` and
`*.stderr.log` for every probe, suite and revert arm.

### 7.1 `atomtypeSevenTest.py`, brought up to the five-leaf reality

The file dated from I-159, when argon had three leaves. It had gone stale in three ways, each now
replaced by a test rather than by an edited comment:

- **its argon closure test covered only `Ar`/`Ar0`/`Ar+`/`Ar++`** — so `Ar0s` and `Ar0e` were
  outside it exactly while this ticket rewired the charge edges *through* them. Now covers all five.
- **no census of the leaf set.** `test_argon_leaf_set_is_exactly_these_five` fails if a sixth leaf
  appears or one of the five goes missing, and checks each is linked to generic `Ar` both ways —
  which per-type membership assertions could not see.
- **"the only two argon signatures in the database" was left implied as complete coverage.** It is
  still factually right (re-verified by grep: `Ar u0 p4 c0` and `Ar u1 p3 c+1` are the only argon
  species in `RMG-database-plasma`), but it means **three of the five leaves — `Ar0s`, `Ar0e`,
  `Ar++` — have no database species at all**, and whatever they do is pinned by unit tests alone.
  `test_argon_leaves_without_a_database_species` states that outright and tells the next ticket to
  move `Ar0e` into `ARGON_DB_SIGNATURES` when a metastable species lands.

  > **`ARGON_DB_SIGNATURES` does not scan a database, and no longer claims to.** It is two literal
  > adjacency lists transcribed by hand from a grep on 2026-09-13, in a file that loads no database
  > at all. So it cannot notice a database that gains a metastable argon species tomorrow; the test
  > above would keep passing. The third round narrowed the constant's comment and both docstrings to
  > that reach: the label is registered, and it is not what either literal types as. Keeping it a
  > literal is the deliberate choice — a scan would need `@pytest.mark.database` and a cloned
  > database at a compatible branch, in a file that is otherwise pure unit test. The cost is that
  > the date is the last moment the completeness claim was checked, and refreshing it is a manual
  > grep. A test that overstates its own reach is the failure mode this campaign keeps recording,
  > so the overstatement was removed rather than the limitation hidden.

It also now pins that `Ar0s` and `Ar0e` produce the *same* sample atom (`Ar u0 p3 c0`) — they differ
only in `single`, the one feature `make_sample_atom` does not act on — which is the shared root
cause of both sitting in `EXPECTED_FAILING_ATOMTYPES` (§5).

---

## 8. What this could not reach

Named explicitly, because most of the value of `Ar0e` is downstream of everything below.

- **No argon species and no plasma family exercises the type.** Both are gated to sibling tickets.
  Everything here is perception and declaration: an atom that can be *built* and *named*. That
  `Ar0e` behaves correctly inside reaction generation — that a template matches it as argon rather
  than as a wildcard, that families produce and consume it, that a metastable survives a
  resonance-structure pass — is **unproven** and cannot be proven from this repository alone.
- **The eleven `resonance.py` sites remain a silent-drop path.** They call
  `update_atomtypes(log_species=False)` and swallow `AtomTypeError` with `pass`. `Ar u2 p3 c0` no
  longer reaches that path, but any *other* untypeable argon still vanishes there without a log
  line. Not fixed, not in scope.
- **The generic-`R` fallback itself is not fixed.** `molecule.py:1587-1593` still assigns the
  wildcard silently when `raise_exception=False`, with logging behind a separate flag.
  `Species.get_resonance_hybrid` (`species.py:765`) is still the one site passing both tolerant
  flags — verified during this ticket, not widened. One atom stopped reaching it; the mechanism is
  untouched.
- **Perception is five times wider than construction — for every argon type, not just `Ar0e`.**
  Because `get_atomtype` ignores `u`, `Ar0e` answers for all five `u` values at `p3 c0`; the
  adjacency list builds one. `Ar0`, `Ar+` and `Ar++` each have the same 5:1 ratio and have had it
  since they were declared (census §2, `test_every_argon_type_perceives_five_u_states`). `Ar0e` did
  add four new perceive-but-not-construct states — `u0`, `u1`, `u3`, `u4` at `p3 c0`, which
  previously raised `AtomTypeError` — and they are reachable mid-recipe: `LOSE_RADICAL` on a
  metastable argon leaves `u1 p3 c0`, which now types as `Ar0e` rather than erroring.

  > **Correction.** The first version of this report called this "the one behaviour I-222 knowingly
  > widens". That was overstated: the pattern is family-wide and pre-existing. `Ar0e` is a fifth
  > instance of it, not the introduction of it.

  Nothing in a declaration can exclude those states — `single`, `lone_pairs` and `charge` are the
  only knobs perception has. Whether any live recipe reaches one was **not** measured; that needs a
  family, which is out of scope.
- **A group spelled `Ar0e` means more than the metastable triplet.** Group adjacency lists get no
  valency-consistency check, so `1 Ar0e ux p3 c0` matches `u1` as readily as `u2`
  (`test_group_spelled_ar0e_matches_more_than_the_metastable_triplet`). `u0` and `u3`+ escape only
  because `update_atomtypes` runs `update_lone_pairs` first and renormalises `p` — arithmetic in a
  different function, not a constraint the label carries. A database group author cannot rely on it.
- **`Ar0e` cannot be built by the sample-molecule machinery** (§5), so any tree-generation or
  group-extension path that relies on `make_sample_molecule` will not produce it. Same limitation
  `Ar0s` has carried since I-218, and for the same reason: the two share a sample atom.
- **Forty-four action slots stay empty where the primitive does produce a type** (§3.3). Those are
  omissions rather than false statements, and were out of scope. Generic `Ar`'s ten self-preserving
  declarations are likewise untouched — they are what make `Ar0e` reactable at all (§3.1), and
  whether that is right is a question about generic `Ar`, not about this leaf.
- **`Ar0e` was chosen and defended, not validated by use.** No group definition, family or library
  spells the label yet, so its ergonomics as an adjacency-list key are demonstrated only by the test
  that writes one.
- **No functional or regression test was run.** Unit and database suites were (§9); functional and
  `test/regression/` were not.

### 8.1 A consequence of this change, referred to the owner

**Not repaired here, and deliberately so: `Plasma_Associative_Ionization_Alkaline_Alkaline` can now
reach metastable argon.** The repair is a database change and that lane is gated; the owner has
taken it. Recorded because it belongs in this change's evidence trail.

> **Status as of the third round: fixed on a sibling database branch, and it leaves a merge-ordering
> constraint that is the owner's to carry.** Both top groups there now read `alkaline u2 px cx`, and
> `alkaline` resolves to Mg and Ca only, so argon cannot match. The constraint: **the database
> narrowing must land before or with this branch**, because this branch alone makes `Ar0e`
> constructible while the *installed* family still has the wildcard top. Nothing in this worktree
> changes, and RMG-database remains untouched by it.

In `RMG-database-plasma/input/kinetics/families/Plasma_Associative_Ionization_Alkaline_Alkaline/groups.py`:

- line 35, top group `A`: `1 *1 R u[2,3,4] px cx` — which `Ar0e` matches on every field;
- line 55, group `B`: `1 *2 R u[2,3,4] px cx`;
- lines 23–28, the recipe: `LOSE_RADICAL *1 2`, `LOSE_RADICAL *2 1`, `FORM_BOND *1 1 *2`,
  `GAIN_CHARGE *1 1`;
- `rules.py` is 8 lines and carries **zero** `entry(` rate rules.

Driven by hand on two metastable argons, that exact recipe produces a well-formed Ar₂⁺ —
`Ar+` bonded to `Ar0s`, net charge +1 — inside a family named for alkaline-earth chemistry:

```
before: ['Ar0e', 'Ar0e']
after:  ['u0 p3 c+1 -> Ar+', 'u1 p3 c0 -> Ar0s']   net charge 1
```

Before `Ar0e`, this was unreachable because metastable argon could not be constructed at all. The
mechanism is §3.1: the family's top group is generic `R`, so the leaf's own `set_actions` never
enter the decision. `test_a_generic_template_matches_and_reacts_metastable_argon` pins the
behaviour, deliberately as behaviour and not as a desideratum — whether that family *should* reach
argon is a database question.

One incidental observation from the same run, not chased: the product molecule kept
`multiplicity 5` from its two `u2` reactants although its atoms end at `u0` and `u1`.
`update_atomtypes` does not recompute multiplicity. Unrelated to argon; noted for whoever meets it.

---

## 9. Wider suites

The brief scoped the count comparison to four files. Both wider suites were run anyway, because an
atom-type addition is exactly the change that can disturb group-tree loading somewhere the four
named files never look. Both are clean.

Both were run in **both rounds**, and the second round's database run is the check that nothing
downstream depended on the four false action edges.

| suite | first round | after the action-graph repair |
|---|---|---|
| unit — `pytest test -m "not functional and not database"` | 3298 passed, 48 skipped (`evidence/full_unit.stdout.log`) | **3317 passed, 48 skipped** (`evidence/full_unit2.stdout.log`) |
| database — `pytest test -m "database"` against `../RMG-database-plasma/input` | 118 passed, 1 xfailed, 2308s (`evidence/db_suite.stdout.log`) | **118 passed, 1 xfailed**, 2595s (`evidence/db_suite2.stdout.log`) |

The third round added no test and changed no executable line, so its job was to show the counts
standing still: **3317 passed, 48 skipped** again after the rebuild
(`evidence/full_unit3.stdout.log`), and 68 / 195 / 69 / 20 collected across the four molecule files
(`evidence/round3_perfile.stdout.log`) — identical to the row below. The database suite was not
re-run in that round: nothing it exercises changed.

Zero failures anywhere. The unit count rises by the 19 new tests. The database count is unchanged at
118 — **no database test changed state when four declared argon edges were corrected**, which is the
evidence that nothing was relying on them. The single `xfail` is
`i134DuplicateElectronsTest::test_one_library_carrying_both_channels_can_be_loaded`, pre-existing.

Risk was also checked statically before the edit: `grep -rln "Ar0\|Ar++" input/kinetics/families/`
over `RMG-database-plasma` returns nothing, so no family spells an argon leaf and none could be
orphaned by dropping `Ar0.increment_charge`.

**One pre-existing collection error, unrelated and not caused here.** Collecting the whole `test`
tree in one pytest run fails before any test executes:

```
ERROR collecting test/rmgpy/rmg/rmgTest.py
import file mismatch: imported module 'rmgTest' has this __file__ attribute:
  test/rmgpy/data/rmgTest.py
which is not the same as the test file we want to collect:
  test/rmgpy/rmg/rmgTest.py
```

Two test files share the basename `rmgTest.py` and neither directory has an `__init__.py`, so
pytest's prepend import mode cannot tell them apart. It reproduces under the repo's own invocation,
`python -m pytest -m "not functional and not database"`, i.e. `make test` verbatim. Both files are
untouched by this branch — `git diff 78f306665 HEAD --` on them is empty, and their last change,
`86652a7d0`, is an ancestor of this branch's base. Clearing `__pycache__` does not help. Both
runs above therefore carry `--ignore=test/rmgpy/rmg/rmgTest.py`; that file was **not** exercised by
either run, which is a gap in this evidence rather than a claim about it. Worth its own ticket: as
things stand, `make test-all` cannot collect this tree in one pass.

The third-round review checked the ignore independently rather than taking it: mainline at
`78f306665` produces the identical error with none of this branch's changes present, and both
colliding files are in the base tree. Recorded here so it is not re-litigated.
