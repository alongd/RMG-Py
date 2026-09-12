# I-218 — Does narrowing `Ar0s` to `single=[1]` remove its test exemption?

Branch `i218-ar0s-single-narrowing`, worktree `/home/alon/Code/RMG-Py-i218-ar0s-single-narrowing`,
base `535c679cf`. Measured 2026-09-12.

## Verdict

**The hypothesis is REFUTED.** Narrowing `Ar0s` from `single=[0, 1]` to `single=[1]` does **not**
remove `'Ar0s'` from `EXPECTED_FAILING_ATOMTYPES`. With the entry removed and the narrowing applied,
`test_make_sample_molecule` still fails, with a *byte-identical* `UnexpectedChargeError` on a
*byte-identical* offending graph. `'Ar0s'` stays in the allowlist.

The ticket also asked whether the two stated consequences are one problem. **They are two.**
Narrowing closes consequence 1 (the `Ar0`/`Ar0s` overlap on the bond-free neutral) and leaves
consequence 2 (the sample builder) exactly as it was. The narrowing is kept on this branch for
consequence 1 alone; it buys nothing toward the allowlist entry.

## 1. Rebuild, and the evidence that it took

The worktree contained **zero** `.so` files at the start — contrary to the brief, which stated that
`rmgpy/molecule/atomtype.cpython-39-x86_64-linux-gnu.so` already existed. Nothing was measured
before a build.

- `python utilities.py check-pydas` → wrote `rmgpy/solver/settings.pxi` with `DEF DASPK = 1`.
- `python setup.py build_ext --inplace` → exit 0, **104** extensions built.
- Every later edit to `atomtype.py` was followed by another `build_ext --inplace`, each of which
  logged `Cythonizing rmgpy/molecule/atomtype.py` and copied a fresh `.so` into `rmgpy/molecule/`.

The freshness assertion is on **values**, not on counts or on the absence of an error:

```
MODULE FILE: /home/alon/Code/RMG-Py-i218-ar0s-single-narrowing/rmgpy/molecule/atomtype.cpython-39-x86_64-linux-gnu.so
IS .so     : True
rmgpy pkg  : /home/alon/Code/RMG-Py-i218-ar0s-single-narrowing/rmgpy/__init__.py
group file : /home/alon/Code/RMG-Py-i218-ar0s-single-narrowing/rmgpy/molecule/group.cpython-39-x86_64-linux-gnu.so
'Ar0s' in ATOMTYPES: True
database.directory: /home/alon/Code/RMG-database-plasma/input
```

This mattered: the environment's `PYTHONPATH` defaults to `/home/alon/Code/RMG-Py/` (the primary
checkout, six trees share it). Every command in this ticket ran with `PYTHONPATH` set to this
worktree, and the printed `__file__` above is the proof it resolved here.

After each edit, `ATOMTYPES['Ar0s'].single` was printed from the *loaded* module before any test —
`[0, 1]` before, `[1]` after. That is the tell that the rebuild reached the binary.

## 2. Correction to the brief: the "233 vs 39" staleness trap does not exist

The brief's Traps section states that `atomtypeTest.py` on a correct build reports **233 passed /
3 skipped**, and that a run reporting **39** is measuring a stale binary. **That is wrong, and my
own baseline refutes it.** On the verified-fresh build above, with `'Ar0s' in ATOMTYPES == True`:

| scope | collected | result |
|---|---|---|
| `atomtypeTest.py` alone | 41 | **39 passed, 2 skipped** |
| `moleculeTest.py` alone | 195 | **194 passed, 1 skipped** |
| both files in one invocation | 236 | **233 passed, 3 skipped** |

39 + 194 = 233. 2 + 1 = 3. 41 + 195 = 236. The campaign's "233 / 3" figure was
`atomtypeTest.py` **plus** `moleculeTest.py`, recorded as if it were a single-file number. The two
figures differ because the **scope** differed, never because a build was stale. A test count cannot
detect a stale `.so` here, and the premise had propagated into three briefs.

The discriminator that does work is the one used above: print the loaded module's `__file__` and
assert on a value the new build must carry.

## 3. Baseline: the failure, reproduced and quoted

With `EXPECTED_FAILING_ATOMTYPES = ["O4b", "S4b"]` and **no other change**, on the unmodified
`single=[0, 1]` declaration:

```
$ python -m pytest test/rmgpy/molecule/atomtypeTest.py::TestAtomType::test_make_sample_molecule -q
E       AssertionError: Couldn't make sample molecules for types Ar0s
E       assert 1 == 0
E        +  where 1 = len(['Ar0s'])
test/rmgpy/molecule/atomtypeTest.py:174: AssertionError

------------------------------ Captured log call -------------------------------
ERROR    root:atomtypeTest.py:172 Couldn't make sample molecule for atomType Ar0s
Traceback (most recent call last):
  File ".../test/rmgpy/molecule/atomtypeTest.py", line 169, in test_make_sample_molecule
    result = group.make_sample_molecule()
  File "rmgpy/molecule/group.py", line 2928, in rmgpy.molecule.group.Group.make_sample_molecule
  File "rmgpy/molecule/group.py", line 3022, in rmgpy.molecule.group.Group.make_sample_molecule
    raise UnexpectedChargeError(graph=new_molecule)
rmgpy.exceptions.UnexpectedChargeError
============================== 1 failed in 1.91s ===============================
```

Instrumenting the same path step by step gives the mechanism and the offending graph:

```
=== make_sample_molecule('1 Ar0s ux') ===
after pick_wildcards, radical_electrons = []
atoms after add_implicit_atoms_from_atomtype = 1
make_sample_atom -> symbol=Ar u=0 p=3 c=0
after update_charge with 0 bonds -> charge = 2
make_sample_molecule RAISED UnexpectedChargeError
offending graph:
1 Ar u0 p3 c+1 {2,S}
2 H  u0 p0 c0 {1,S}
```

## 4. Why narrowing cannot fix that — measured, then read off the source

Re-running the identical probe against the rebuilt `single=[1]` binary produced output **identical
in every line but the declaration**:

```
$ diff baseline_probe.txt narrowed_probe.txt
2c2
< Ar0s single=[0, 1] lone_pairs=[3] charge=[0]
---
> Ar0s single=[1] lone_pairs=[3] charge=[0]
```

Same `atoms after add_implicit_atoms_from_atomtype = 1`, same `u=0 p=3 c=0`, same `charge = 2`, same
`UnexpectedChargeError`, same `1 Ar u0 p3 c+1 {2,S}` + H graph. And `pytest` on the narrowed build
with the entry removed still reports `1 failed`.

The reason, in the source:

1. `Group.add_implicit_atoms_from_atomtype` (`group.py:2452`) adds implicit partners for
   `o_double`, `s_double`, `r_double`, `all_double`, `triple`, `quadruple` and lone pairs. It
   **never adds a single-bonded partner**. Narrowing `single` to `[1]` therefore does not cause a
   bonded sample to be built — the atom stays bond-free.
2. `GroupAtom.make_sample_atom` (`group.py:818`) takes `radical_electrons[0]` if the group has any
   and the `Atom` default (`0`) otherwise. A `ux` group carries none, so the sample is `u0`. A
   neutral argon at `p3` needs `u1` to balance; nothing in the builder can choose it.
3. The saturation loop (`group.py:2978-2984`) clips added hydrogens to
   `max(atomtype.single) - single_present`. `max([0, 1])` and `max([1])` are **both 1**. Either way
   one H is added, leaving `c+1`, and `Ar0s` is not on the `positive_charged` allowlist at
   `group.py:3002`, so `UnexpectedChargeError` is raised.

So the sample-builder failure is structurally independent of the `single` list. Removing the
allowlist entry needs the builder to pick a radical count that balances the charge — a change in
`make_sample_atom`, which this ticket explicitly did not commission. Per the brief, that is where I
stopped.

## 5. What narrowing *does* fix — consequence 1

Both builds were driven through every bond-free neutral argon the perceiver can be handed:

| adjacency list | `single=[0, 1]` (baseline) | `single=[1]` (narrowed) |
|---|---|---|
| `1 Ar u0 p3 c0` | `InvalidAdjacencyListError: ... for atom Ar (Ar0s) ...` | `InvalidAdjacencyListError: ... for atom Ar (Ar) ...` |
| `1 Ar u1 p3 c0` | `InvalidAdjacencyListError: ... for atom Ar (Ar0s) ...` | `InvalidAdjacencyListError: ... for atom Ar (Ar) ...` |
| `1 Ar u2 p3 c0` | **`Ar0s`** | **`AtomTypeError: Unable to determine atom type`** |
| `1 Ar u0 p4 c0` | `Ar0` | `Ar0` |

Two things to read here. First, the error text on the top two rows names the atom type the
perceiver reached before the adjacency-list consistency check rejected it: `(Ar0s)` on the baseline,
plain `(Ar)` on the narrowed build. That is the brief's own trap made visible — `get_atomtype()`
ignores radical electrons, so `Ar0s` was claiming all three `u` states at zero bonds, and only the
later valency check was refusing two of them. Second, `1 Ar u2 p3 c0` — the bond-free triplet, the
one state that *did* balance — typed as `Ar0s` on the baseline and is refused outright after
narrowing. That is the overlap closing: after the narrowing, the bond-free neutral argon belongs to
`Ar0` alone, and it belongs there by declaration rather than by lone-pair count and `specific`-list
ordering.

**Why closing it is safe, on stronger evidence than "no test failed."** The manager checked, and I
confirmed independently, that the `u2 p3 c0` state the narrowing makes untypeable is *uninhabited*:
`RMG-database-plasma/input/` contains exactly two argon species, `Ar u0 p4 c0` (37 occurrences,
types `Ar0`) and `Ar u1 p3 c+1` (4 occurrences, types `Ar+`). There is no bond-free argon at `p3`
anywhere in the database — no argon metastable species exists in it at all; PlasmaAir's "metastable
quenching" header refers to the N/O metastables. Narrowing therefore orphans no real species.
(The *absence* of an Ar\* metastable is a database question the manager is filing separately; it is
not chased here.)

## 6. Controls

All three ran on both builds and are unchanged by the narrowing.

**A — `Ar2+` still builds and still types `['Ar0s', 'Ar+']`.** The reason `Ar0s` exists.

```
Ar2+ atom types = ['Ar0s', 'Ar+']
Ar2+ net charge = 1
multiplicity    = 2
CONTROL A PASS   (both builds)
```

**B — a bare argon atom still types `Ar0`, unambiguously.** Checked through `update_atomtypes()`
*and* through `get_atomtype()` directly, and through group matching:

```
bare argon type          = Ar0
get_atomtype(bare Ar)    = Ar0
'1 Ar u0 p4 c0' matches '1 Ar0 ux'  : True
'1 Ar u0 p4 c0' matches '1 Ar0s ux' : False
CONTROL B PASS   (both builds)
```

**C — `Ar+` / `Ar++` resolution unchanged.**

```
Ar u1 p3 c+1 -> Ar+
Ar u0 p3 c+2 -> Ar++
CONTROL C PASS   (both builds)
```

## 7. Test counts, before and after

Each file in its own scope, compared against its own collection count. Allowlist entry present
(`["O4b", "S4b", "Ar0s"]`) in all rows.

| suite | collected | before (`single=[0,1]`) | after (`single=[1]`) |
|---|---|---|---|
| `test/rmgpy/molecule/atomtypeTest.py` | 41 | 39 passed, 2 skipped | **39 passed, 2 skipped** |
| `test/rmgpy/molecule/moleculeTest.py` | 195 | 194 passed, 1 skipped | **194 passed, 1 skipped** |
| `test/rmgpy/molecule/groupTest.py` | 69 | 69 passed | **69 passed** |
| `test/rmgpy/molecule/atomtypeSevenTest.py` | 11 | 11 passed | **11 passed** |

313 passed, 3 skipped across the four suites, identical before and after. The narrowing is inert on
this surface.

The reproduction rows, where the allowlist entry is removed, are the point of the whole ticket —
and they too are unchanged by the narrowing:

| build | scope run | result |
|---|---|---|
| `single=[0, 1]` | `atomtypeTest.py::TestAtomType::test_make_sample_molecule` | 1 failed |
| `single=[1]` | `atomtypeTest.py::TestAtomType::test_make_sample_molecule` | 1 failed |
| `single=[1]` | `atomtypeTest.py` (full file) | **1 failed, 38 passed, 2 skipped** |

The full file with the entry removed was run on the narrowed build only; on the baseline build the
single test node was run. Both failures are the same `UnexpectedChargeError` quoted in §3.

## 8. What this may and may not be claimed to show

**Measured: direct atom-type perception, and a narrow slice of group matching.**

- Direct perception — `get_atomtype()` and `Molecule.update_atomtypes()` over every bond-free
  neutral argon state, plus `Ar+`, `Ar++` and both halves of `Ar2+`.
- Group matching — `AtomType.is_specific_case_of()` for `Ar0s` against all six of its generic
  labels (`R`, `R!H`, `R!H!Val7`, `Rx`, `Rx!H`, `Ar`; all `True`, unchanged by the narrowing), and
  `Molecule.is_subgraph_isomorphic()` of `Ar2+` and of a bare argon against the `Ar0s`, `Ar0` and
  `R` group patterns (unchanged).

**NOT measured: the database group/family matching surface.** The brief's warning stands and is not
answered by this work. Adding `Ar0s` to the broad generic lists widened group matching through
`equivalent()` / `is_specific_case_of()`; `groupTest.py` and the molecule suites exercise the
*machinery* of group matching on hand-built groups, not the RMG-database group trees, kinetics
family templates, or thermo group estimation. No reaction family was driven, no database was
loaded, no model was generated. A narrowing that is inert on 313 unit tests can still move which
template a real family node matches.

**This is therefore not a merge-ready result**, and nothing here should be read as verifying the
branch. It answers one question — narrowing does not lift the allowlist entry — and closes one
declaration/name mismatch. Clearing it for merge needs an adversarial round over the family- and
thermo-group matching surface.

Also unverified by construction: the `u` states `Ar0s` admits were measured through the adjacency
list consistency check, which is what actually rejects them, not through the atom type — as the
brief required. The atom type itself still ignores radical electrons.

## 9. What changed on this branch

```
rmgpy/molecule/atomtype.py   single=[0,1] -> single=[1] on Ar0s, plus a comment block recording
                             why it is [1] and why that does NOT lift the allowlist entry
docs/i218-ar0s-single-narrowing/  this report, stdout.log, stderr.log
docs/contracts/i218-ar0s-single-narrowing.md
```

`EXPECTED_FAILING_ATOMTYPES` is **unchanged** — it was edited only for the reproduction in §3 and
restored. No action edge was added or removed on any atom type. `Ar0`, `Ar+`, `Ar++`, `Mg0s` and
`Ca0s` are untouched, as is `atomtypeSevenTest.py`. `rmgpy/solver/settings.pxi` is build output and
is not committed. Nothing was pushed and nothing was merged.
