# Round 107 — the fields survive the loader and die in the transforms

Base `i221-saturation-no-atomtype@5a821898d`. Gates: `rmgpy/data/kinetics/family.py`,
`rmgpy/data/kinetics/library.py`, and the test tree.

Both findings confirmed. The HIGH's census is wider than the review states: the class has
**nine** members, not three, and **two of the ones in gates were missed** — one of which
is where round 105's own three fields were being dropped all over again.

---

## 1. HIGH — the transforms enumerated by hand, and the census is nine sites

### 1.1 What was measured

`docs/i221-saturation/probes/round107_transform_enumeration_probe.py`, run at
`5a821898d` (`logs/round107-probe-red-at-5a821898d.log`). Every field below was set on
the reaction before the transform and read back after it; the control
`the fixture holds every field before the transform` proves the left-hand column.

| transform | field | before → after |
|---|---|---|
| `TemplateReaction.__reduce__` | `elementary_high_p` | `True → False` |
| | `allow_pdep_route` | `True → False` |
| | `allow_max_rate_violation` | `True → False` |
| | `rank` | `7 → None` |
| | `comment` | `'a comment worth keeping' → ''` |
| | `label` | `'Lip + CH3 <=> CH3Li' → ''` |
| | `network_kinetics` | `3.0 → None` |
| `TemplateReaction.copy()` | the same three flags | `True → False` |
| | `rank` | `7 → None` |
| | `comment` | **`'…' → None`** |
| | `network_kinetics` | `3.0 → None` |
| | `labeled_atoms` | `['products', 'reactants'] → []` (unset) |
| `LibraryReaction.__reduce__` | `rank` | `7 → None` |
| | `comment` | `'…' → ''` |
| | `label` | `'…' → ''` |
| | `is_forward` | **`True → False`** |
| `LibraryReaction.copy()` | `allow_max_rate_violation`, `rank`, `is_forward` | reset |
| | `entry` | `'an entry' → None` |
| | `library`, `family` | **gone — the attribute does not exist** |
| | *the class* | **`LibraryReaction → Reaction`** |

Three things the review did not state:

- **`copy()` leaves `comment` as `None`**, not `''`. `Reaction` declares it
  `cdef public str`; `TemplateReaction.copy()` built the object with `__new__` and never
  assigned the field, so the copy holds a value the class says it cannot hold. Every
  `reaction.comment.split(...)` downstream — including `authoring_family`, which reads
  provenance out of exactly that string — raises on it rather than returning nothing.
- **`LibraryReaction.__reduce__` was dropping round 105's three fields.** Round 105 closed
  `is_forward`, `rank` and `comment` disappearing in `as_library_reaction`. They
  disappear again the first time the result is pickled. Same defect, same cause, a
  different hand-written list, one transform over.
- **`LibraryReaction` had no `copy()` at all.** It inherited `Reaction.copy`, which does
  `Reaction.__new__(Reaction)`, so a copy of a library reaction was a *base reaction*:
  `library`, `family` and `entry` had nowhere to land. `entry` is the carrier the
  quarantine gate reads authorship from, which round 95 added to the `TemplateReaction`
  pair for precisely this reason and nobody carried across to the twin. It is reached:
  `rmgpy/tools/isotopes.py:394` copies whatever `Reaction` it is handed, and the core
  model it walks is full of library reactions.

### 1.2 The census — nine sites, four in gates

| site | in gates | state it dropped at `5a821898d` |
|---|---|---|
| `reaction.py:322` `Reaction.__reduce__` | no | `allow_max_rate_violation`, `is_forward` |
| `reaction.py:2017` `Reaction.copy` | no | the same two; and returns a base `Reaction` for every subclass that does not override it |
| `family.py` `TemplateReaction.__reduce__` | **yes** | 7 fields, table above — **repaired** |
| `family.py` `TemplateReaction.copy` | **yes** | 7 fields — **repaired** |
| `library.py` `LibraryReaction.__reduce__` | **yes** | 4 fields — **repaired** |
| `library.py` `LibraryReaction.copy` | **yes** | everything the subclass adds — **added** |
| `library.py` `_carry_entry_fields` | **yes** | nothing; rounds 102/105 made it derived |
| `depository.py:83` `DepositoryReaction.__reduce__` | no | the three flags, `rank`, `comment`, `label`, `network_kinetics` |
| `rmg/pdep.py:90` `PDepReaction.__reduce__` | no | the three flags, `rank`, `comment` |
| `rmg/model.py:79` `as_library_reaction` | no (this round) | nothing; round 105 made it derived |

Ten rows, nine distinct sites plus the already-derived conversion. The review named three.
The two it missed inside the gates are the `LibraryReaction` pair. `depository.py` and
`pdep.py` are the same defect in files this round cannot touch: **named, not fixed.**

### 1.3 The repair, and why this shape

One derived definition, in `family.py`, shared by all five in-gates sites:

```python
def state_fields(source):
    return REACTION_STATE_FIELDS.union(getattr(source, '__dict__', {}))
```

Two halves, because a `Reaction` subclass has two kinds of field. `Reaction` is a Cython
extension type and its fields are descriptors on the class, which is what
`REACTION_STATE_FIELDS` has discovered since round 102. `TemplateReaction`,
`LibraryReaction` and `DepositoryReaction` are **plain Python classes**, so the fields
they add (`family`, `template`, `estimator`, `reverse`, `entry`, `labeled_atoms`;
`library`) are ordinary instance attributes, invisible to a `dir()` of the class and
visible in `vars()` of the instance. Measured: `vars()` gives exactly those six on a
`TemplateReaction`, those three on a `LibraryReaction`, and **nothing** on a base
`Reaction`, which has no instance dictionary at all — so for the loader and the
conversion, which pass the default universe, nothing changes.

`__reduce__` now returns `(cls, (), state)` — **no constructor arguments at all**. That is
the point rather than a detail. A `__reduce__` that rebuilds through a constructor can only
carry what the constructor accepts, and `TemplateReaction.__init__` has no parameter for
any of the three flags; keeping a constructor-argument list beside a complete state dict
would leave a hand-written list in the file that no longer matters, which is the thing that
rots. With no arguments there is one definition and nothing to drift.

`__setstate__` is required, not cosmetic. The default unpickling path is
`inst.__dict__.update(state)`, and every `Reaction` field is a data descriptor on the
class, which takes precedence over the instance dictionary — so the values would land
somewhere nothing ever reads. The explicit `__setstate__` goes through `setattr`, and
through `_CARRIED_THROUGH_STORAGE` for `degeneracy`, so round 105's rule (carrying a value
is not the same operation as changing it) is inherited by the transforms rather than
restated in them.

`copy()` keeps its deep half explicit and derives its shallow half. The exclusion table
names each deepened field and why; the control `no transform restates the rate while
carrying the degeneracy` holds in **both** directions, so routing `degeneracy` through the
generic helper did not reintroduce round 105's rescale.

### 1.4 What changed that a reader should know

The pickle format changed shape. A pickle written **before** this tip still loads here: it
carries no third element, so `__setstate__` is never called and the old constructor path
runs. A pickle written **here** needs the new `__setstate__`. RMG pickles reactions within
a run — multiprocessing, `deepcopy` — and writes seed mechanisms as `.py` libraries, so
there is no on-disk format to migrate. Stated in the probe as a note rather than worked
around.

---

## 2. MEDIUM — the carry claimed to fail loudly and did not

`carry_reaction_state` caught `AttributeError` from both the source read and the target
write and continued. Its docstring said a field could "only be left behind by someone
writing down why". Both statements cannot be true.

Measured at `5a821898d`, two ways:

- take `protons` out of the loader's exclusion policy — claiming it carried — and the
  carry **returns normally**. `protons` is read-only (control: assigning it raises), so
  every such field was a silent drop wearing the label of a carry. That is what `protons`
  itself was for the whole of round 102.
- hand the carry a source object holding none of the state at all, and it **completes in
  silence**, returning a reaction it wrote nothing to. Every `getattr` raised, every one
  was swallowed.

It raises `ReactionStateNotCarried` now, from both sides, naming the field and the class.

**Deliberately not an `AttributeError` subclass.** Both underlying failures *are*
`AttributeError`s, and these call sites sit inside loaders that catch `AttributeError`
liberally — an exception one of those could swallow would restore exactly the silence the
raise removes. There is a test for that, because it is the kind of thing a later
"tidy-up" undoes.

Raising is safe because the partition makes it unreachable, and the partition is pinned:
`test_every_field_called_carried_can_be_carried_between_the_real_classes` runs the real
carry over every policy × every **(source class, target class) pair used in production** —
`Reaction → TemplateReaction`, `Reaction → LibraryReaction`, `TemplateReaction →
LibraryReaction`, and each class to itself. Round 105's version probed a bare `Reaction`
against itself, which is not a pair any carry runs on; it has been deleted rather than
kept alongside, with a comment at its old site saying why.

---

## 3. The evidence note — structural vs behavioural red states

The note is accepted and acted on in two places.

**Relabelled in place.** Round 105's four structural reds now carry the label in their own
docstrings, naming the base they are structural at and the behavioural test that covers the
same ground. The convention is written into the test module's docstring, so it applies to
the next round without being restated. A findings document goes stale; a docstring beside
the test does not.

**This round's split, measured** (`logs/round107-tests-red-at-5a821898d.log`): 17 of the 19
new tests fail at `5a821898d`, 164 pass, nothing skips — the same 181 the tip reports, so
the comparison is field-for-field.

- **8 behavioural.** `test_the_three_flags_and_the_carrier_survive` in all four
  shape × transform combinations, failing with the dropped field by name; the three MEDIUM
  refusals, failing with `DID NOT RAISE` or "carried in silence"; and the
  docstring-agreement test. None of these imports a name the repair adds — that is
  deliberate, and it is why they can run at the base at all.
- **9 structural.** The partition, discovery and census tests, failing with `ImportError`
  for `_NOT_REPRODUCED`, `state_fields` or `_TEMPLATE_NOT_COPIED_BY_REFERENCE`. They pin
  the shape of the repair. Each names its behavioural counterpart in its docstring.

The ratio is worse than round 105's (17/4) and that is honest rather than accidental: this
round's repair *is* largely a shape, and the behaviour it fixes is covered by 8 tests plus
7 probe findings rather than by 17. Two tests are **green at base by construction** —
`test_the_copy_does_not_restate_the_rate[template|library]`, guards on the new carry path,
which the probe's control measures in both directions.

Probe: **7 of 7 findings reproduce at `5a821898d`, 0 of 7 at the tip; 6 of 6 controls hold
at both.**

---

## 4. What I could not reach

- **Whether any of this has ever bitten a real run.** The four transforms are proven wrong
  by construction and by measurement, and `isotopes.py:394` is a live caller of `copy()`,
  but I did not run an isotope job or a pressure-dependent deck to watch a dropped
  `elementary_high_p` change a mechanism. The defect is certain; its historical blast
  radius is not.
- **`depository.py` and `rmg/pdep.py`.** Both drop the same fields, both are out of gates,
  both are measured in the probe's census note. `PDepReaction`'s `__reduce__` matters most
  of the three remaining, because pressure-dependent networks are pickled across processes
  by design.
- **`Reaction.copy` returning a base `Reaction`.** Fixed for `LibraryReaction` by giving it
  an override; `DepositoryReaction` still has none, so a copy of one is still a base
  reaction with `depository`, `family` and `entry` gone. `rmgpy/reaction.py` is out of
  gates, and the general repair belongs there.
- **Whether `SurfaceArrhenius` and `SurfaceChargeTransfer` should exist.** They are
  `cdef public` slots named after cimported types, assigned by nothing in the codebase and
  always `None`. Round 105's exclusion reason called them "cimported types"; measured, they
  are never-assigned slots. The reason strings are corrected in `family.py` and
  `library.py`; **the same two strings in `rmgpy/rmg/model.py:104` still say "a cimported
  type"** and are out of this round's gates, so three policies now carry two wordings for
  one fact. Removing the slots from `reaction.pxd` altogether is also out of gates and is
  probably the right answer.
- **Performance.** Building a state dict per pickle is not measured. It is dominated by the
  `Species` deep copies already in the same tuple, but that is an argument, not a number.

## 5. Named, not fixed

- `rmgpy/data/kinetics/depository.py:83` and `rmgpy/rmg/pdep.py:90` — §1.2.
- `rmgpy/reaction.py:322` and `:2017` — the base pair, dropping `allow_max_rate_violation`
  and `is_forward`, and returning a base `Reaction` from `copy()`.
- `rmgpy/tools/isotopes.py:458` — round 105's HIGH, still live: `degeneracy` assigned three
  lines after `rxn.kinetics = new_kinetics`, so an isotopologue's rate is still rescaled.
- `rmgpy/data/kinetics/database.py:493` — drops `electrons`. Recommended as the next HIGH
  for the fourth round running.
- `report.md` §19.2 (database repo) names engine commits only through `9b89a937c` and is
  now five rounds stale. Update it as the last act before any merge, not before.
