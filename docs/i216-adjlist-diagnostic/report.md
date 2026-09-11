# I-216 — adjlist valency error destroys its own diagnostic

Worktree: `/home/alon/Code/RMG-Py-i216-adjlist-diagnostic`, branch `i216-adjlist-diagnostic`.

## Environment resolution

- `rmgpy.settings['database.directory']` resolved to `/home/alon/Code/RMG-database/input`
  (default fallback #5 in the priority chain — no `rmgrc` exists in this worktree or `~/.rmg/`).
  The directory exists on disk. No test in this session actually loaded the thermo/kinetics
  database (adjacency-list parsing and `ConsistencyChecker` do not touch it), so this resolution
  is reported for completeness per the brief, not because it was load-bearing.
- Import provenance, confirmed via `inspect.getsourcefile`/`inspect.getsource`:
  - `rmgpy.molecule.adjlist` -> this worktree's `rmgpy/molecule/adjlist.py` (pure Python, not
    cythonized -- no `adjlist*.so` exists or is produced by `make build`). Edits take effect
    immediately, no rebuild needed.
  - `rmgpy.molecule.atomtype.get_atomtype` -> this worktree's compiled
    `rmgpy/molecule/atomtype.cpython-39-x86_64-linux-gnu.so`.
- The worktree initially had zero `.so` files under `rmgpy/` (`ModuleNotFoundError: No module
  named 'rmgpy.rmgobject'` on first import -- see stderr.log). Fixed by running `make build`
  (incremental in-place Cython build; safe, scoped to this checkout, does not touch the shared
  `rmg_env`). No `pip install -e .` or `make unsafe-install-shared-env` was run.

## Step 1 -- reproducing the masking (before touching anything)

Claim under test: `get_atomtype(atom, atom.edges).label`, called inline inside the
`InvalidAdjacencyListError(...)` message construction in `check_partial_charge`, can raise
`AtomTypeError` and this exception propagates in place of the valency diagnostic.

Driven through the real `Molecule().from_adjacency_list(...)` path with:
```
1 Ar u0 p4 c0 {2,S}
2 Ar u0 p4 c0 {1,S}
```
An Ar atom given any bond cannot match any of Ar's specific atom types (`Ar0`/`Ar+`/`Ar++`, all
defined with zero bond-forming actions), so `get_atomtype` raises `AtomTypeError` for it -- and
this atom also fails the valency check (bonded Ar with 4 lone pairs is not charge-balanced).
Exact traceback obtained (quoted verbatim from `docs/i216-adjlist-diagnostic/stdout.log`):

```
Traceback (most recent call last):
  File "<string>", line 8, in <module>
  File "rmgpy/molecule/molecule.py", line 1942, in rmgpy.molecule.molecule.Molecule.from_adjacency_list
    def from_adjacency_list(self, adjlist, saturate_h=False, raise_atomtype_exception=True,
  File "rmgpy/molecule/molecule.py", line 1951, in rmgpy.molecule.molecule.Molecule.from_adjacency_list
    self.vertices, self.multiplicity, self.metal, self.facet = from_adjacency_list(adjlist, group=False, saturate_h=saturate_h,
  File "/home/alon/Code/RMG-Py-i216-adjlist-diagnostic/rmgpy/molecule/adjlist.py", line 873, in from_adjacency_list
    ConsistencyChecker.check_partial_charge(atom)
  File "/home/alon/Code/RMG-Py-i216-adjlist-diagnostic/rmgpy/molecule/adjlist.py", line 109, in check_partial_charge
    type=get_atomtype(atom, atom.edges).label,
  File "rmgpy/molecule/atomtype.py", line 1041, in rmgpy.molecule.atomtype.get_atomtype
    def get_atomtype(atom, bonds):
  File "rmgpy/molecule/atomtype.py", line 1069, in rmgpy.molecule.atomtype.get_atomtype
    raise AtomTypeError(
rmgpy.exceptions.AtomTypeError: Unable to determine atom type for atom Ar, which has 1 single bonds, 0 double bonds (0 to O, 0 to S, 0 others), 0 triple bonds, 0 quadruple bonds, 0 benzene bonds, 4 lone pairs, and +0 charge.
```

**Result: CONFIRMED, no contradiction.** The exception that escapes `from_adjacency_list` is
`rmgpy.exceptions.AtomTypeError`, raised inside the `.format(...)` call at `adjlist.py:109`,
*instead of* the intended `InvalidAdjacencyListError` at line 112 -- the valency diagnostic is
never constructed, let alone raised. This matches the ticket's claim exactly.

## Step 2 -- the fix

Minimal, idiomatic diff to `rmgpy/molecule/adjlist.py`:

```diff
-from rmgpy.exceptions import InvalidAdjacencyListError
+from rmgpy.exceptions import AtomTypeError, InvalidAdjacencyListError
@@ check_partial_charge
         if not (-0.301 < atom.charge - theoretical < 0.301):
             # It should be 0, but -0.1 is caused by a Hydrogen bond
+            try:
+                atom_type_label = get_atomtype(atom, atom.edges).label
+            except AtomTypeError:
+                # The atom type can't be perceived (often *because* of the bad valency
+                # being reported here). Fall back to the element symbol so the valency
+                # error itself isn't masked by an unrelated AtomTypeError.
+                atom_type_label = atom.symbol
             raise InvalidAdjacencyListError(
                 'Invalid valency for atom {symbol} ({type}) with {radicals} unpaired electrons, '
                 '{lone_pairs} pairs of electrons, {charge} charge, and bonds [{bonds}].'.format(
                     symbol=atom.symbol,
-                    type=get_atomtype(atom, atom.edges).label,
+                    type=atom_type_label,
                     radicals=atom.radical_electrons,
                     lone_pairs=atom.lone_pairs,
                     charge=atom.charge,
                     bonds=','.join([str(bond.order) for bond in atom.bonds.values()])
                 )
             )
```

Only `AtomTypeError` is caught (the exact exception the reproduction showed escaping) -- no bare
`except`, no change to any valency rule, no new/modified atom type.

## Step 3 -- the two controls (real pytest tests)

Added to `test/rmgpy/molecule/adjlistTest.py`, class `TestConsistencyChecker`:
- `test_check_partial_charge_bad_valency_untypeable_atom` -- Ar dimer as above.
- `test_check_partial_charge_bad_valency_typeable_atom` -- control: a carbon with 4 single bonds
  to H (features match ordinary `Cs`) but declared `u2` (2 unpaired electrons), which
  `get_atomtype`'s feature matching ignores (it does not look at `radical_electrons`), so the atom
  types cleanly as `Cs` while still failing the valency check.

RED confirmation (pre-fix, source reverted to HEAD via a scoped `git stash` isolating only
`adjlist.py`, tests already present):
```
FAILED test/rmgpy/molecule/adjlistTest.py::TestConsistencyChecker::test_check_partial_charge_bad_valency_untypeable_atom
1 failed, 1 passed, 31 deselected in 2.02s
```
(the control test passed even pre-fix -- expected, since it never enters the code path that
raises `AtomTypeError`.)

GREEN confirmation (post-fix):
```
test_check_partial_charge_bad_valency_untypeable_atom PASSED
test_check_partial_charge_bad_valency_typeable_atom PASSED
2 passed, 31 deselected in 1.55s
```

Exact messages captured post-fix:
```
Untypeable atom:  'Invalid valency for atom Ar (Ar) with 0 unpaired electrons, 4 pairs of electrons, 0 charge, and bonds [1.0].'
Typeable control: 'Invalid valency for atom C (Cs) with 2 unpaired electrons, 0 pairs of electrons, 0 charge, and bonds [1.0,1.0,1.0,1.0].'
```

Byte-identity of the control message pre-fix vs post-fix: the control case never raises
`AtomTypeError` (get_atomtype succeeds, returns `Cs`) in either version of the code, so the
`try/except` is never entered -- the two versions execute the identical `.format(...)` call with
identical inputs. This is not just an assumption: the control test (asserting the exact string
`Invalid valency for atom C (Cs)`) passed unmodified against both the pre-fix and post-fix source,
which is the operational form of "byte-identical" for this case.

## Step 4 -- sibling sweep (report only, not fixed)

Searched all `raise InvalidAdjacencyListError(...)` / `.format(...)` sites in `adjlist.py`
(~35 locations). All of them format plain data already in hand -- strings (`atom.symbol`, `line`,
the raw `adjlist` text), ints (`aid`, `atom1`, `atom2`, `multiplicity`, `n_rad`) -- none of which can
raise. The only other function call anywhere near error-message construction is
`get_element(atom.number, isotope)` at line ~805, but that is a plain attribute assignment
(`atom.element = get_element(...)`), not part of building an error message, so it is out of scope
for this defect shape.

**Finding: no other instance of the same defect shape (a call that can raise, invoked inline
inside error-message construction) exists in `adjlist.py`.** Nothing else was fixed, per
non-goals.

## Step 5 -- full-suite counts

Before (unfixed `adjlist.py`, restored via targeted `git stash`; new tests already present):
```
791 passed, 14 skipped, 1 deselected in 7.84s     (full suite, minus the untypeable test)
1 failed, 1 passed, 31 deselected in 1.74s         (the two new tests alone: untypeable test FAILED as expected, control PASSED)
```
Combined "before" total: 806 collected items -> 791 passed + 14 skipped + 1 failed (the untypeable
test, which is expected to fail pre-fix by design) = 806. No other failures were observed in the
before run -- no pre-existing failures unrelated to this ticket were found.

After (fixed `adjlist.py`, stash re-applied and dropped):
```
792 passed, 14 skipped in 7.92s
```
792 = 791 + 1 (the previously-failing test now passes); 14 skipped unchanged; 0 failed. Exit
code 0.

## What this may and may not be claimed to show

**May claim:** the masking described in the ticket is real and was reproduced exactly, for at
least one concrete, minimal case (a bonded Ar atom failing both atom-type perception and the
valency check). The fix restores the intended `InvalidAdjacencyListError` valency diagnostic in
that case, naming the atom by element symbol when its type cannot be perceived, while leaving the
message for an ordinary typeable atom with a bad valency unchanged. The full `test/rmgpy/molecule/`
suite has no regressions (792 passed, 14 skipped, same skip count before and after).

**May NOT claim:** that this is the only element/bonding combination capable of triggering the
same masking (only Ar was used as the concrete witness; the argument for other untypeable-atom
cases is structural, not exhaustively tested), or that no other latent
call-that-can-raise-inside-error-message-construction defect exists anywhere else in the RMG-Py
codebase -- the sweep in Step 4 was scoped to `adjlist.py` only, as instructed.

## Contradictions with the brief

None found. Every empirical check (masking reproduction, control construction, RED-before-fix,
GREEN-after-fix, byte-identity of the control, before/after suite counts) matched what the ticket
predicted.
