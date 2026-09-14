# A check that could not fail, and the false positive that repairing it exposed

Branch `databaseTest-cd-atom-type-repair`, off `plasma@311818121`. Two commits, both in
`test/database/databaseTest.py` and one new test file. Nothing outside the test tree is touched.

This branch is intended to travel with **I-230**, which installs a different new check into the
same file. Three independent changes to `databaseTest.py` arriving as three separate collisions is
the avoidable cost.

## What was wrong

Both Cd/CO/CS/Cdd checks — `kinetics_check_cd_atom_type` (the 140 kinetics families) and
`general_check_cd_atom_type` (called from `test_thermo`, `test_solvation`, `test_statmech` and
`test_transport`, covering 19 group databases) — decided "is this a double bond?" like this:

```python
num_of_d_bonds = sum([1 if x.order[0] == "D" and len(x.order) == 1 else 0 for x in atom.bonds.values()])
```

A `Group` bond stores `order` as a list of *numbers*. `order[0] == "D"` was never true,
`num_of_d_bonds` was always `0`, and **neither check could fail for anything it was ever given**.
The branch immediately below it in the same function has always used the numeric test; they now
agree.

## What repairing it found

Exactly one complaint across the whole database — and it was the check, not the data. The
solvation node `O2d-Cdd` is a carbon declared `Cdd`, which carries *two* double bonds, with one of
them drawn. **A group is a subgraph**, so `num_of_d_bonds` counts the double bonds that are drawn
and is only a lower bound. Reading it as the total, the check concluded "exactly one, and it is to
an O, so this must be `CO`" — and `CO` and `Cdd` are mutually exclusive atom types, neither a
specific case of the other, so the demand was not merely unmet but unsatisfiable.

The second commit narrows that inference to fire only when the atom is not already declared `Cdd`.
The `num_of_d_bonds == 2 ⇒ Cdd` branch is untouched and stays sound for a subgraph, because a
carbon cannot carry more than two.

## Which findings a test pins, and which are only prose

| # | finding | pinned by | if it regresses |
|---|---|---|---|
| 1 | the repaired predicate actually fires — a `Cd` double-bonded to `O` and not typed `CO` is caught, at both call sites | `cdAtomTypeTest.py::TestCdAtomTypeCheckIsLive` (4 tests) | **red** |
| 2 | the reason it was dead: a `Group` bond order is numeric, never the letter `"D"` | `TestCdAtomTypeCheckCannotRegress` (2 tests). The first is a characterization of `Group` and fails *naming both checks* if that ever changes | **red** |
| 3 | the narrowing is sound and bought only the `Cdd` case — a `Cd`-typed atom with one drawn double bond to O is still caught, two drawn bonds still demand `Cdd`, and `CO`/`Cdd` really are mutually exclusive | `TestCdAtomTypeInferenceIsSoundForSubgraphs` (6 tests, including `test_the_allowance_did_not_switch_the_check_off`) | **red** |
| 4 | **no entry in the live database fails the repaired checks** — 140 kinetics families and 19 group databases | pinned by `databaseTest.py` itself: `test_kinetics`, `test_thermo`, `test_solvation`, `test_statmech`, `test_transport` all run the repaired predicate over the real data. **6 passed in 342.27s** (`logs/06-databaseTest-narrowed.stdout.log`) | **red** |
| 5 | the *size* of that zero — 17333 `Group` entries scanned, 15775 atoms typed Cd/CO/CS/Cdd, 12851 of them double-bonded | **prose + probe only.** `probes/cd_blast_radius.py`, run by hand (`logs/07-blast-radius-final.stdout.log`). Nothing in the suite re-runs it. | silent |
| 6 | the tests were confirmed RED before the repair went in — 2 failed, 4 passed against an unrepaired copy taken from `311818121` | **prose + log only** (`logs/02-red-arm.stdout.log`). This is a claim about the history of the work, and no test can pin it. | n/a |

Row 5 is the one to know about. The *conclusion* it supports — nothing in the database fails the
repaired checks — is row 4 and is pinned live by the suite; only the reach numbers behind it are
prose. That distinction matters because a zero from a check that has never run is otherwise
indistinguishable from a check that still does not run, which is why the probe reports reach at
all, and why it refuses to print anything until a negative control confirms the as-shipped
predicate is silent exactly where the repaired one complains.

**No pin in this file asserts an enumeration order.** The one positional expression,
`group.atoms[0]` in `test_the_old_string_comparison_would_have_been_dead`, indexes into a group
built from a literal adjacency list inside the same test, so the position is authored rather than
incidental.

## The probe's own asymmetry, and why

The two arms of `probes/cd_blast_radius.py` are deliberately unlike each other. The **repaired**
arm imports `databaseTest.py` by path and calls the real check functions, so it cannot drift away
from what the suite runs. The **as-shipped** arm is re-implemented inside the probe, because that
code no longer exists in the file. A probe that re-implements *both* arms is a mirror with nothing
to tie it back to the implementation, and this campaign has shipped that mistake before.

## Runs

| run | arm | result |
|---|---|---|
| `test/database/cdAtomTypeTest.py` | unrepaired copy of `databaseTest.py` from `311818121` | 2 failed, 4 passed — the RED arm |
| `test/database/cdAtomTypeTest.py` | repair only, before the narrowing | 6 passed |
| `test/database/databaseTest.py` | repair only, before the narrowing | **1 failed**, 5 passed — the `O2d-Cdd` false positive |
| `test/database/cdAtomTypeTest.py` | repair + narrowing (HEAD) | 12 passed |
| `test/database/databaseTest.py` | repair + narrowing (HEAD) | 6 passed in 342.27s |

Both streams are captured for every run under `logs/`.
