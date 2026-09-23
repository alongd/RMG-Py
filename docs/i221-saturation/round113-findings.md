# Round 113: four partial closures, and a green that went stale

Base: `i221-saturation-no-atomtype@fc60e5ba4`. Commits: `a30f61d1f` (item 1), `520f6bda3` (item 2),
`040b54279` (item 3), `d66b9a7a3` (item 4). All are local. Nothing was pushed, merged or rebased.
Database: `/home/alon/Code/RMG-database-plasma/input` at `0d9c5bc2462847c631ed0e977ed2a1dd1c3b4648`,
read-only. Interpreter: `/home/alon/anaconda3/envs/rmg_env/bin/python`.

Each item was reproduced red at `fc60e5ba4` before its fix. Logs are in `logs/round113-*`. Both
streams were captured for each run.

## 1. HIGH 2: the reverse's sides are now rebuilt

`make_new_reaction` (`rmgpy/rmg/model.py`) rewrote `forward.reverse.pairs` over the model species.
It left `forward.reverse.reactants/products` as the objects the family built, which are Molecules
from `add_reverse_attribute`. So no reverse pair member occurred on its own side. The reverse sides
are now the forward's model species, swapped, in fresh lists. The reverse's own template order is
not kept, because its Molecules carry no link to the model species and the reverse pairs already
followed the forward's order.

The old test compared the rebuilt pairs with their own reversal. It now asserts
`pair_occurrences(reverse.pairs, reverse.reactants, reverse.products) == [(1, 1), (0, 0), (2, 0)]`
and per-slot identity. A second test builds the reverse from Molecules in a different order, as
production does.

- Red (`logs/round113-item1-red.log`, `-k Model`): **2 failed, 3 passed**. The failure is the
  brief's `ValueError: ... is not a reactant of this reaction`.
- Green (`logs/round113-item1-green.log`, the whole file): **19 passed**.

**Not reached:** the end-to-end probe's reactions are not own-reverse, so the probe does not
exercise this path. The coverage is the unit test.

## 2. HIGH 5: any unanswered loaded family makes the enumeration incomplete

`_enumerate_quarantines` recorded an unanswered loaded family only when its cached `quarantine`
was set. The dropped case was a label the disk pass then skipped as already seen. The premise
behind the old condition was that a None attribute is unanswered only when no families directory
exists. That premise was false: `resolve_quarantine` answers a loaded family from its attribute in
exactly that case. Every unanswered resolution is now recorded, whatever the attribute holds.

There are two new tests. The first is the reviewer's reproduction: `_family_directory` returns
`None`, `listdir` returns `['F']`, and `F.quarantine = None`. The second puts a loaded family with
a None attribute behind a denied directory.

The permission tests used `chmod 000` and skipped when the process could still read, which is the
case under root. They now inject `PermissionError` by patching `os.stat`, `os.lstat`, `os.listdir`
and `os.scandir` for the one directory, through `_deny_access`. `os.open` is left alone:
`_read_manifest` requires `os.open in os.supports_dir_fd`, and a wrapper fails that check.

- Red (`logs/round113-item2-red.log`): **2 failed, 11 passed**. The failures are `([], [])` and
  `DID NOT RAISE`.
- Green (`logs/round113-item2-green.log`, the whole file): **261 passed, 1 failed**. The one
  failure is item 4, which was still open at that point.
- Run as uid 0 under `unshare -r`: the `fc60e5ba4` test file gave **3 passed, 8 skipped**
  (`logs/round113-item2-root-skips-at-r112.log`). The current file gives **13 passed, 0
  skipped** (`logs/round113-item2-root-now.log`).

## 3. MEDIUM: a subclass defined after install is refused

`install_complete_reducers` registered a refusal in `ForkingPickler`'s table for each subclass that
existed when it ran. That table is keyed on the exact type, so a later subclass fell through to its
inherited reducer and was rebuilt as the parent class. The refusal now resolves by MRO at reduce
time, through one `reducer_override` (`_refuse_by_mro`) that `_CompletePickler` and `ForkingPickler`
share. Verdicts are cached per class in a `WeakKeyDictionary`. The per-subclass table entries are
gone. `install` refuses to overwrite a `reducer_override` that someone else installed.

The cost is one cached lookup per non-builtin object on the multiprocessing path. Round 112 had
chosen a table entry to avoid that cost, and that choice is what left the hole.

The new test defines the subclass after install, then sends it through `ForkingPickler.dumps` and a
real `Pipe`.

- Red (`logs/round113-item3-red.log`): **1 failed, 11 passed**. The failure is `DID NOT RAISE
  PicklingError`.
- Green (`logs/round113-item3-green.log`, subclass, transport and census tests): **43 passed**.

## 4. The `registered` fixture read the live database

The fixture replaced the loaded registry but left `settings['database.directory']` live. The
on-disk half of enumeration therefore read `Plasma_Electron_Impact_Ionization/quarantine.py`, an
Arrhenius manifest, which produced two warnings in a control that must stay silent. The fixture
now points the setting at an empty families directory of its own, and the control asserts that.
The production warning is untouched.

**Was round 112 false-green? No.** I checked by running the `fc60e5ba4` test against a `git
archive` of the database at `96f2afa4a`: **1 passed**. Against `0d9c5bc24`: **1 failed**
(`logs/round113-item4-premise.log`). The manifest arrived at DB merge `ceb49fcb7` (2026-09-23).
Round 112's count was true against the database it ran on. The fixture had never been hermetic,
and the database drift exposed that.

In the sweep, `TestTheOtherTransportGetsTheSameReducers` and `TestASpeciesOnBothSidesKeepsItsSide`
load real families through the live setting, but they carried no `database` marker. They are now
marked, like their siblings. Under `-m database` the file gives **25 passed**
(`logs/round113-item4-database-marked.log`).

- Red (`logs/round113-item4-red.log`): **1 failed**.
- Green (`logs/round113-item4-green.log`, the whole file): **263 passed**.

## Verification

Focused suite, the brief's command, against the database above
(`logs/round113-focused-suite.log`): **510 passed, 4 skipped (WIP markers), 49 deselected, 0
failed.** How that reconciles with round 112's 519:

| step | count |
|---|---|
| round 112's reported pass count (518 passed + item 4's failing test) | 519 |
| new tests (items 1, 2, 2, 3) | +4 |
| tests moved behind the `database` marker (13 in the two classes) | −13 |
| **this run** | **510** |

End-to-end probe `probes/round112_end_to_end_probe.py` (`logs/round113-e2e.log`): **49 ok, 0
failures**.

The five pre-existing `i134DuplicateElectronsTest` failures are outside the focused paths and were
not touched.
