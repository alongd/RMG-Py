# I-222 — An atom type for bond-free neutral argon (metastable Ar\*)

Branch `i222-metastable-argon-atomtype`, worktree
`/home/alon/Code/RMG-Py-i222-metastable-argon-atomtype`, base `78f306665`.
Nothing pushed, nothing merged.

**Outcome: the type is added.** `Ar u2 p3 c0` — metastable argon — now builds as a `Molecule` and
as a `Species` through the ordinary path and types as **`Ar0e`**. The tolerant path no longer
degrades it to the generic wildcard `R`.

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

---

## 3. `set_actions`

`Ar0e` declares **no** action edges — all ten lists empty, as `Ar0s` does. This is not a default; it
is what the measurement forces. Applying each single action to a concrete `u2 p3 c0` argon atom
(`evidence/actions.stdout.log`):

| action                | resulting state | types as   |
|-----------------------|-----------------|------------|
| `GAIN_PAIR`           | `u2 p4 c-2`     | *no type*  |
| `LOSE_PAIR`           | `u2 p2 c+2`     | *no type*  |
| `LOSE_CHARGE`         | `u2 p3 c-1`     | *no type*  |
| `GAIN_CHARGE`         | `u2 p3 c+1`     | **`Ar+`**  |
| `FORM_BOND` (single)  | `u2 p3 c0`, 1 bond | **`Ar0s`** |
| `GAIN_RADICAL`        | `u3 p3 c0`      | `Ar0e`     |
| `LOSE_RADICAL`        | `u1 p3 c0`      | `Ar0e`     |

So every action that lands anywhere at all lands on a **sibling**. Declaring
`increment_charge=['Ar+']` needs `Ar+` to declare `decrement_charge=['Ar0e']` back; declaring
`form_bond=['Ar0s']` needs `Ar0s` to declare `break_bond=['Ar0e']`. Both inverse entries live on
atom types this ticket may not modify, and one side alone is precisely the one-way edge
`TestActionGraphClosure` refuses. This is the same reason `Ar0s` declares nothing, recorded there
for I-218.

The only edges declarable without touching a sibling are the radical **self**-edges, which would
close trivially. They are omitted for consistency: `Ar0`, `Ar+` and `Ar++` each map to themselves
under `GAIN_RADICAL` too, and each declares `increment_radical=[]`. The whole argon family omits
them.

Because no edge is added, the action graph is unchanged and all three closure tests stay green,
including `test_argon_and_alkali_families_close_both_ways`, which fails on any argon entry landing
in the pre-existing allowlist.

### Two pre-existing declarations the probe exposed, and did not touch

Both are out of scope (`Ar0`, `Ar+` are named non-goals) and neither is affected by this change.
Recorded because the probe measured them and the next argon ticket will meet them:

- `ATOMTYPES['Ar0'].set_actions(..., increment_charge=['Ar+'], ...)` — but `GAIN_CHARGE` on a real
  `Ar0` atom (`u0 p4 c0`) gives `u0 p4 c+1`, which types as **nothing at all**, not `Ar+`.
- `ATOMTYPES['Ar+'].set_actions(..., decrement_charge=['Ar0'], ...)` — but `LOSE_CHARGE` on a real
  `Ar+` atom (`u1 p3 c+1`) gives `u1 p3 c0`, which types as **`Ar0e`** now (and as nothing before
  this change), not `Ar0`. `Ar0` is at `p4`; a bare charge action does not move lone pairs.
- Likewise `Ar0.decrement_lone_pair = ['Ar+']`, where `LOSE_PAIR` on `Ar0` actually gives
  `u0 p3 c+2` = `Ar++`.

These are label symmetries that satisfy the closure test while naming states the actions do not
produce. Closure checks that the graph is symmetric, not that its edges are true.

---

## 4. What changed, by file

`rmgpy/molecule/atomtype.py`
- new `ATOMTYPES['Ar0e']` declaration and its all-empty `set_actions`, each with the reasoning above
  as comments;
- `'Ar0e'` added to `ATOMTYPES['Ar'].specific` and to the `specific` lists of `R`, `R!H`,
  `R!H!Val7`, `Rx`, `Rx!H`.
- `Ar0`, `Ar0s`, `Ar+`, `Ar++` are **unmodified**; `git diff` shows no line of any of them changed.

`test/rmgpy/molecule/atomtypeTest.py`
- `EXPECTED_FAILING_ATOMTYPES` gains `"Ar0e"` (see §5). `"Ar0s"` stays.
- the two tripwires rewritten in place and renamed, per §0:
  `test_bond_free_triplet_argon_has_no_atom_type` → `test_bond_free_triplet_argon_types_as_ar0e_not_ar0s`,
  `test_untypeable_argon_degrades_to_generic_R_when_typing_is_tolerant` →
  `test_metastable_argon_no_longer_degrades_to_generic_R_when_typing_is_tolerant`;
- new class `TestMetastableArgonAtomType`, 12 tests.

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

`rmgpy/molecule/atomtype.py` was restored to `git show 78f306665:...` — tests left in place — and
the extension rebuilt. The loaded module was checked by value, not by mtime:

```
LOADED: .../rmgpy/molecule/atomtype.cpython-39-x86_64-linux-gnu.so
Ar0e in ATOMTYPES: False
```

Result (`evidence/revert_red.stdout.log`): **13 failed, 43 passed, 2 skipped**. All 11 `Ar0e`
assertions in the new class fail, and so do both rewritten tripwires:

```
FAILED TestArgonSingleBondNarrowing::test_bond_free_triplet_argon_types_as_ar0e_not_ar0s
FAILED TestArgonSingleBondNarrowing::test_metastable_argon_no_longer_degrades_to_generic_R_when_typing_is_tolerant
FAILED TestMetastableArgonAtomType::test_declaration_is_the_bond_free_p3_neutral
FAILED TestMetastableArgonAtomType::test_metastable_argon_builds_as_a_molecule
FAILED TestMetastableArgonAtomType::test_metastable_argon_builds_as_a_species
FAILED TestMetastableArgonAtomType::test_tolerant_typing_yields_ar0e_and_not_the_wildcard
FAILED TestMetastableArgonAtomType::test_no_two_argon_types_can_match_one_atom
FAILED TestMetastableArgonAtomType::test_each_realizable_argon_atom_types_specifically
FAILED TestMetastableArgonAtomType::test_only_u2_is_constructible_at_p3_c0
FAILED TestMetastableArgonAtomType::test_perception_admits_unrealizable_u_states
FAILED TestMetastableArgonAtomType::test_ar0e_declares_no_action_edges
FAILED TestMetastableArgonAtomType::test_ar0e_is_linked_into_the_type_hierarchy_both_ways
FAILED TestMetastableArgonAtomType::test_ar0e_is_usable_as_a_group_adjacency_list_label
```

The twelfth new test, `test_untypeable_argon_still_degrades_to_the_wildcard`, **passed** against the
reverted tree — as it must. It is the preserved coverage of the generic-`R` fallback, and its whole
value is that it holds on both sides of the change.

`atomtype.py` was then restored and rebuilt; `atomtypeTest.py` returned to 56 passed / 2 skipped.

---

## 7. Suite counts, per file against its own collected total

| file | collected (base) | result (base) | collected (now) | result (now) |
|---|---|---|---|---|
| `atomtypeTest.py`      | 46  | — | **58** | 56 passed, 2 skipped |
| `moleculeTest.py`      | 195 | — | **195** | 194 passed, 1 skipped |
| `groupTest.py`         | 69  | — | **69** | 69 passed |
| `atomtypeSevenTest.py` | 11  | — | **11** | 11 passed |
| four together          | 321 | 318 passed, 3 skipped | **333** | **330 passed, 3 skipped** |

`atomtypeTest.py` grows by exactly 12 — the new class. No other file's collected total moves, and no
test anywhere goes from passing to failing. The two skips in `atomtypeTest.py` are the pre-existing
`@pytest.mark.skip(reason="WIP")` sample-molecule tests; the third skip is
`moleculeTest.py::test_count_internal_rotors_dimethyl_acetylene`, also pre-existing.

Full unit suite (`pytest test -m "not functional and not database"`): see §9.

Both streams were captured for every measurement; `evidence/` holds the matching `*.stdout.log` and
`*.stderr.log` for the baseline probe, the action probe, the baseline suite, the post-change suite
and the revert run.

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
- **A widening this change knowingly introduces.** Because perception ignores `u`, `Ar0e` also
  answers for `u1 p3 c0` and `u3 p3 c0` — states no molecule can hold, since the adjacency list
  refuses them. Before this change they raised `AtomTypeError`. They are reachable mid-recipe:
  `LOSE_CHARGE` on `Ar+` leaves `u1 p3 c0`, which now types as `Ar0e` instead of erroring. Nothing
  in the declaration can exclude them — `single`, `lone_pairs` and `charge` are the only knobs
  perception has. `test_perception_admits_unrealizable_u_states` pins this as a hazard so that a
  later cost has a test to read rather than a surprise to debug. Whether any live recipe reaches it
  was **not** measured, because that needs a family, which is out of scope.
- **`Ar0e` cannot be built by the sample-molecule machinery** (§5), so any tree-generation or
  group-extension path that relies on `make_sample_molecule` will not produce it. Same limitation
  `Ar0s` has carried since I-218.
- **The three pre-existing false action edges in §3 were measured but not fixed**, per the non-goal.
  One of them, `Ar+.decrement_charge = ['Ar0']`, now names a *different* wrong type than it did
  before this change — the action produces `Ar0e`, where before it produced nothing typeable. The
  declaration was wrong both ways; this change moved which way.
- **No database, functional or regression test was run.** The type is new, so nothing in
  RMG-database references it yet; but that is an argument, not a measurement.
- **`Ar0e` was chosen and defended, not validated by use.** No group definition, family or library
  spells the label yet, so its ergonomics as an adjacency-list key are demonstrated only by the test
  that writes one.

---

## 9. Wider suites

The brief scoped the count comparison to four files. Both wider suites were run anyway, because an
atom-type addition is exactly the change that can disturb group-tree loading somewhere the four
named files never look. Both are clean.

**Unit suite** — `pytest test -m "not functional and not database"`
(`evidence/full_unit.stdout.log`):

```
3298 passed, 48 skipped, 153 deselected, 63 warnings in 153.77s
```

**Database suite** — `pytest test -m "database"`, against `../RMG-database-plasma/input`
(`evidence/db_suite.stdout.log`):

```
118 passed, 3380 deselected, 1 xfailed, 2 warnings in 2308.78s
```

Zero failures in either. The single `xfail` is
`i134DuplicateElectronsTest::test_one_library_carrying_both_channels_can_be_loaded`, pre-existing.

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
