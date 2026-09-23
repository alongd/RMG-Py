# Round 108 — copying the reaction component by component severed its internal references

Base `i221-saturation-no-atomtype@56156ac9a`. Gates: `rmgpy/data/kinetics/family.py`,
`rmgpy/data/kinetics/library.py`, and the test tree.

The HIGH is confirmed and closed, and both halves of it — the defect and the tests'
structural blindness to it. Two corrections to the brief, both measured: **`pickle` never
had this defect**, so one of the five acceptance criteria is a property the base already
satisfies; and **the suggested repair does not work for half the problem**, because
`Molecule.__deepcopy__` discards the memo a shared-memo `deepcopy` would depend on.

---

## 1. The finding, measured

`docs/i221-saturation/probes/round108_copy_identity_probe.py` at `56156ac9a`
(`logs/round108-probe-red-at-56156ac9a.log`) — **4 of 8 findings reproduce, 4 of 4
controls hold**:

| transform | `pairs` member owned by the copy | labelled atom inside a copied species | relabel observable |
|---|---|---|---|
| `copy()`, template | **no** | **no** | **no** |
| `copy()`, library | **no** | *(no `labeled_atoms`)* | — |
| `pickle`, template | yes | yes | yes |
| `pickle`, library | yes | — | — |

The defect is real and its consequences are the two the review names:

- `Species.__eq__` is **identity** (`species.py:205`, `return self is other`), so a
  `pairs` entry that is not one of the copy's own species makes
  `reactants.index(pair[0])` raise `ValueError` — loudly, but far from the copy.
- `labeled_atoms` is read back at `family.py:2593` to relabel the reaction's structures
  immediately before `pairs` and `template` are regenerated. A detached atom makes that a
  no-op: the labels go onto atoms nothing else holds, the structures stay unlabelled, and
  the regenerated pairs and template are derived from them. **Silent**, and at the one
  point where being wrong is invisible.

### What round 107 did and did not cause

I take the criticism, with one correction that the measurement supports.

- **`pairs`**: the severing is **older than round 107**. `other.pairs = deepcopy(self.pairs)`
  sat beside separately-copied reactants in `TemplateReaction.copy` before I touched it,
  and it still sits in `Reaction.copy` (`reaction.py:2017`) today. Round 107 carried it
  across into `LibraryReaction.copy`, which is new exposure for that class but not a new
  defect.
- **`labeled_atoms`**: **this one is mine.** Before round 107 `copy()` did not set the
  field at all, so reading it raised `AttributeError`. I gave it a detached deep copy. A
  loud absence became a silent wrong value, which is strictly worse, and it is the second
  time on this branch that widening a carry introduced a quieter failure than the one it
  replaced.

The through-line the review states is right, and it has a specific shape worth naming: the
repairs have each been correct about *what* to carry and naive about *what carrying means*.
Round 105 — carrying a value through a setter that transforms it. Round 108 — carrying a
reference by cloning what it points at. Both times the field arrived and the meaning did
not.

---

## 2. Correction — `pickle` never had this defect

The Verifier asks for the same four assertions after a `pickle` round trip. They hold at
`56156ac9a` already, in every combination, and the probe records that as four findings
that do **not** reproduce rather than four that the repair fixed.

The reason is structural: `__reduce__` returns one state dict, and pickle's memo spans the
whole of it, so two references to one species are two references to one restored species by
construction. That was true before round 107 as well, when the same fields travelled as one
constructor-argument tuple.

The four pickle tests are kept — a property worth holding is worth pinning — but they are
labelled in their own docstrings as green at the base, not as evidence. Presenting them as
part of the repair would be claiming a fix for something that was never broken.

---

## 3. Correction — a shared-memo `deepcopy` cannot do this, and the reason is upstream

The review offers two shapes and I tried the second one first, because it is the better
one: *"use a single `deepcopy` with a shared memo so the cross-references are preserved by
construction."*

Measured, through one memo:

| | shared with the copied species? |
|---|---|
| a `pairs` **species** | **yes** |
| a `labeled_atoms` **atom** | **no** |

`rmgpy/molecule/molecule.py:1064`:

```python
def __deepcopy__(self, memo):
    return self.copy(deep=True)
```

The override accepts the memo and discards it. Every molecule is therefore rebuilt with
fresh atoms no matter what the caller has already copied and recorded, so **no caller-side
memo discipline can preserve atom identity through a `deepcopy`**. Species-level sharing
survives only because `Species` has no such override. `rmgpy/molecule/` is out of gates:
named, not changed.

`pickle` keeps its own memo, which `__deepcopy__` never sees. So the repair reproduces the
deepened half through one pickle round trip — the same mechanism `__reduce__` already
defines, which makes `copy()` and the pickle path one implementation rather than two that
have to agree. It is also **3× faster** than the `deepcopy` it replaces (0.130 ms against
0.406 ms on a four-species reaction; the four `Species.copy(deep=True)` calls it replaces
cost 0.116 ms), so the correctness is not bought with time.

There is a negative control on that choice:
`test_a_deepcopy_cannot_do_this_and_the_reason_is_upstream` asserts that `deepcopy` still
*fails* to share. If upstream ever honours the memo, it goes red and says the simpler
mechanism has become available — so the workaround cannot outlive its reason.

---

## 4. The repair

Both `copy()` bodies are one call to `copy_reaction(self, <policy>)`. The helper:

1. carries the shallow half with `carry_reaction_state` — `entry` stays the shared
   database object, deliberately, and a control pins that;
2. deepens the rest in **one** round trip, over the set computed as
   `not_copied_by_reference - _NOT_REPRODUCED`.

The two halves are complementary *in code*, not by agreement, so a field added to either
table lands in exactly one of them. And the two `copy()` methods are one implementation
rather than two that agree — pinned by `test_the_two_copies_are_one_implementation`, which
also refuses a `deepcopy` reappearing in either body.

One test of round 107's went red on the refactor:
`test_the_sites_that_enumerate_reaction_state_share_one_definition` required each site to
name one of the three helpers, and a site that delegates to `copy_reaction` names none of
them. That was the checker being too literal — delegation to the shared helper *is* sharing
the definition — so `copy_reaction` was added to the list it accepts. Recorded here rather
than fixed quietly, because "the guard refused my construct" is exactly the moment to say
which of the two was wrong.

---

## 5. The tests, and why they could not see this

The review's second half is right and worth stating precisely. Round 107's fixtures used
`pairs = [("Lip", "CH3Li")]` and `labeled_atoms = {"reactants": {"*1": "Lip"}}` — strings.
A string cannot point into a copied structure, so **no assertion over those fixtures could
fail on a severed reference**, however thorough. Worse, they compared by *value*, and a
detached clone has every right value: the probe carries that as a control
(`value comparison passes whether or not the references survive`) which holds at both tips.

Both fixtures now hold real `Species` and a real `Atom` taken out of a real `Molecule`, and
every assertion in the new class is `is`. The class also carries an anti-vacuity test
asserting the fixture holds the invariant *before* anything is copied — because the string
version did not, and that is the failure mode being repaired.

**Red state** (`logs/round108-tests-red-at-56156ac9a.log`): 14 new tests, of which **5 fail
at the base and 9 pass**.

- **4 behavioural** — the `copy` arms of the pairs, labelled-atom and relabelling tests,
  each failing with the severed reference by name.
- **1 structural but base-runnable** — `test_the_two_copies_are_one_implementation` fails
  on the base's actual source (two bodies, each with its own `deepcopy` calls) rather than
  on a missing import. It pins the repair's shape; the four above carry the defect.
- **9 green at the base by construction** — the four `pickle` arms (§2), the anti-vacuity
  fixture check, the two shallow-half controls, the `Molecule.copy` ordering pin, and the
  `deepcopy` negative control. Each says so in its docstring, per round 107's convention.

---

## 6. What I could not reach

- **Whether a real run has ever consumed a severed copy.** `isotopes.py:394` calls
  `copy()` on whatever reaction it is handed and is the live caller, but I did not run an
  isotope job and watch a template regenerate from unlabelled structures. The mechanism is
  proven; the historical damage is not.
- **`deepcopy(reaction)` is still severed, and it is a live path.** `family.py:3854` and
  `database.py:755` both `deepcopy` a reaction. Because `deepcopy` recurses into
  `Molecule.__deepcopy__`, labelled atoms are severed there exactly as they were in
  `copy()` — and I cannot fix it from inside these gates without either changing
  `molecule.py` or giving `Reaction` a `__deepcopy__`, which belongs in `reaction.py`.
  **This is my recommendation for the next round**, and it is a one-line fix in
  `rmgpy/reaction.py` (`__deepcopy__ = lambda self, memo: self.copy()`) whose correctness
  now rests on a `copy()` that is already right.
- **`Reaction.copy` itself** (`reaction.py:2017`) still severs `pairs` for every subclass
  without an override — `DepositoryReaction` and `PDepReaction` today.
- **Cost in a real generation loop.** The 3× speedup is a microbenchmark on one reaction,
  not a profile of a run.
- **Whether `labeled_atoms` values are ever lists in practice.** `family.py:2594` handles
  `isinstance(atom, list)`, so the tests handle it, but I found no path that produces one
  and did not establish whether that branch is live or vestigial.

## 7. Named, not fixed

- `rmgpy/molecule/molecule.py:1064` — `Molecule.__deepcopy__` discards the memo. The root
  cause of §3, and the reason `copy()` cannot use `deepcopy`.
- `rmgpy/reaction.py:2017` — `Reaction.copy` severs `pairs` and returns a base `Reaction`.
- `rmgpy/data/kinetics/family.py:3854`, `rmgpy/data/kinetics/database.py:755` —
  `deepcopy(reaction)`, severed labelled atoms, live.
- `rmgpy/data/kinetics/depository.py:83`, `rmgpy/rmg/pdep.py:90` — round 107's
  hand-enumerated `__reduce__` pair, unchanged.
- `rmgpy/tools/isotopes.py:458` — round 105's HIGH, still live.
- `rmgpy/data/kinetics/database.py:493` — drops `electrons`; fifth round named.
- `rmgpy/rmg/model.py:104` — round 107's two stale exclusion reasons, out of gates there.
- `report.md` §19.2 (database repo) — six rounds stale. Update before any merge, not
  before.
