# i213 — Does one atom type for a bonded neutral argon make Ar₂⁺ buildable?

**Verdict: yes, and the boundary is sharper than expected.** One atom type — `Ar0s`, seven added
lines in `rmgpy/molecule/atomtype.py` — makes `Ar₂⁺` build, atom-type, survive `update()`, round-trip
through its adjacency list, and balance in a trial reaction. It does **not** make the species usable:
the SMILES path corrupts its charge, and group additivity refuses it three separate ways (loudly —
no silent number). The change is a measuring instrument; nothing here argues for keeping it.

Branch `i213-bonded-argon-atomtype`, worktree
`/home/alon/Code/RMG-Py-i213-bonded-argon-atomtype`, base `546d1c727`. Nothing pushed, nothing merged.

## Provenance

| | |
|---|---|
| `cwd` | `/home/alon/Code/RMG-Py-i213-bonded-argon-atomtype` |
| `python` | `/home/alon/anaconda3/envs/rmg_env/bin/python` |
| `rmgpy.__file__` | `…/RMG-Py-i213-bonded-argon-atomtype/rmgpy/__init__.py` |
| `atomtype.__file__` | `…/rmgpy/molecule/atomtype.cpython-39-x86_64-linux-gnu.so` |
| `atomtype compiled` | `True` |
| **`database.directory`** | **`/home/alon/Code/RMG-database-plasma/input`** (exists; realpath identical) |

The worktree had **no `.so` files and no `rmgrc`** on arrival. `rmgrc` was created from
`rmgrc.template` (git-ignored, database pinned at the plasma checkout beside this tree). The tree was
built with `make build` **before** the "before" measurement and **again after** the source edit, so
neither column was taken against a stale binary. The staleness trap was closed by measurement rather
than by trust: the probe script's provenance block asks the *loaded binary* which argon labels it
holds and prints them (`['Ar', 'Ar+', 'Ar++', 'Ar0']` before, `[… , 'Ar0s']` after).

Reproduce with:

```bash
conda activate rmg_env
cd /home/alon/Code/RMG-Py-i213-bonded-argon-atomtype
python docs/i213-bonded-argon-atomtype/probe.py \
  > >(tee -a docs/i213-bonded-argon-atomtype/stdout.log) \
  2> >(tee -a docs/i213-bonded-argon-atomtype/stderr.log >&2)
```

Logs in this directory: `before-*.log` (probe on the unmodified branch), `stdout.log`/`stderr.log`
(probe after the change), `tests-before-*.log`, `tests-after-noallowlist-*.log`, `tests-after-*.log`,
`widening-*.log`. The `make build` transcripts were discarded — 118 kB of gcc warnings that carry
nothing the provenance block above does not state more directly.

## 1. The original failure, reproduced

Adjacency list under test — the localised-charge form of Ar₂⁺ (ground state ²Σu⁺, a doublet cation):

```
1 Ar u0 p3 c+1 {2,S}
2 Ar u1 p3 c0 {1,S}
```

On the unmodified branch, `Molecule().from_adjacency_list` of that string raises:

```
Traceback (most recent call last):
  File ".../docs/i213-bonded-argon-atomtype/probe.py", line 85, in <lambda>
    lambda: Molecule().from_adjacency_list(AR2_PLUS_ADJLIST),
  File "rmgpy/molecule/molecule.py", line 1953, in rmgpy.molecule.molecule.Molecule.from_adjacency_list
    self.update_atomtypes(raise_exception=raise_atomtype_exception)
  File "rmgpy/molecule/molecule.py", line 1592, in rmgpy.molecule.molecule.Molecule.update_atomtypes
    raise
  File "rmgpy/molecule/molecule.py", line 1587, in rmgpy.molecule.molecule.Molecule.update_atomtypes
    atom.atomtype = get_atomtype(atom, atom.edges)
  File "rmgpy/molecule/atomtype.py", line 1069, in rmgpy.molecule.atomtype.get_atomtype
    raise AtomTypeError(
rmgpy.exceptions.AtomTypeError: Unable to determine atom type for atom Ar., which has 1 single
bonds, 0 double bonds (0 to O, 0 to S, 0 others), 0 triple bonds, 0 quadruple bonds, 0 benzene
bonds, 3 lone pairs, and +0 charge.
```

That is the ticket's error verbatim (the brief drops the radical dot in `atom Ar.`; nothing else
differs). The traceback frames matter: the failure is in `get_atomtype` reached from
`update_atomtypes`, **not** in the adjacency-list consistency checker — see §7, where a different
argon case turns out to be masked by that checker.

### The ticket's correction is confirmed — the limit was never structural

`update_atomtypes` stops at the first atom it cannot resolve, so the two centres were typed
individually (`raise_atomtype_exception=False`, then `get_atomtype` per atom):

| environment | before the change |
|---|---|
| bare `Ar u0 p4 c0` | `Ar0` |
| bare `Ar u1 p3 c+1` | `Ar+` |
| **bonded `Ar+` centre of Ar₂⁺** | **`Ar+` — types today, no change needed** |
| **bonded neutral `Ar` centre of Ar₂⁺** | **`AtomTypeError`** |

So the engine was one atom type short of the species, not structurally incapable of it. The record's
"structural limit" claim is refuted by measurement, exactly as the ticket said.

## 2. The change

Seven lines added to `rmgpy/molecule/atomtype.py`, plus `'Ar0s'` inserted in the `Ar` generic's
`specific=[...]` and in the five element rosters that already name `'Ar0','Ar+','Ar++'`
(`Rx`, `Rx!H`, `R`, `R!H`, `R!H!Val7` — lines 263, 279, 313, 327, 344; all five found by grep, none
guessed):

```python
ATOMTYPES['Ar0s'] = AtomType('Ar0s', generic=['R', 'R!H', 'R!H!Val7', 'Rx', 'Rx!H', 'Ar'], specific=[],
                            single=[0,1], all_double=[0], r_double=[0], o_double=[0], s_double=[0],
                            triple=[0], quadruple=[0], benzene=[0], lone_pairs=[3], charge=[0])
```

### Why `lone_pairs=[3]` and not `Ar0`'s `[4]`

The reconciliation the brief asked for, done as arithmetic rather than assumed. Neutral argon brings
**8** valence electrons. One covalent bond consumes **one** of them (the partner supplies the other),
leaving **7** non-bonding electrons on the atom: **3 lone pairs plus 1 unpaired electron**. Not 4
lone pairs — 4 lone pairs plus a bond would be 9 electrons. So `lone_pairs=[3]`, and the neutral
bonded argon is necessarily a **radical**. That is not an artefact of the model: it is why Ar₂⁺ is a
doublet, and it is the fact that §5 turns into a test failure.

The same arithmetic checks the two neighbours already in the table. `Ar+` has 7 valence electrons:
1 (bond) + 6 (`p3`) = 7 with `u0`, or 6 (`p3`) + 1 radical = 7 unbonded — hence `single=[0,1]`,
`lone_pairs=[3]`, `charge=[1]`, which is what the table declares.

`single=[0,1]`: one bond is the case that matters; `0` is retained because the sibling entries all
admit their bond-free case and because the type must be reachable when a recipe breaks the bond.
`charge=[0]` by definition — this is the *neutral* centre.

**Name.** `Ar0s` follows `Mg0s` / `Ca0s`: the digit is the charge, the `s` says the bonds are single.
It does not follow the halogen convention (`Cl1s`, `F1s`, where the digit counts bonds), because
consistency inside the argon block — where `Ar0`/`Ar+`/`Ar++` are charge-labelled — matters more.

### `set_actions`: deliberately inert, and what that costs

```python
ATOMTYPES['Ar0s'].set_actions(increment_bond=[], decrement_bond=[], form_bond=[], break_bond=[],
                              increment_radical=[], decrement_radical=[], increment_lone_pair=[],
                              decrement_lone_pair=[], increment_charge=[], decrement_charge=[])
```

Two edges would be chemically right and are **not** declared:

- `Ar0s --increment_charge--> Ar+` (the neutral half of Ar₂⁺ losing its electron), and
- `Ar0 --form_bond--> Ar0s`.

Each needs its inverse declared back on `Ar+` / `Ar0`, i.e. an edit to an atom type this ticket may
not touch. Declaring one side alone opens a one-way edge, which `TestActionGraphClosure` in
`test/rmgpy/molecule/atomtypeTest.py` refuses by design. **Consequence to carry forward: no reaction
recipe can move an atom into or out of `Ar0s`.** A species carrying it can be written down and read
back; it cannot be *reached* by a family recipe. If this probe is ever promoted to a change worth
keeping, those two inverse pairs are the first thing to add.

## 3. What now builds

All from `stdout.log`, after `make build`:

| measurement | result |
|---|---|
| `Molecule().from_adjacency_list(Ar₂⁺)` | builds — `<Molecule "[Ar][Ar+]">` |
| atom types | `['Ar0s', 'Ar+']` |
| `mol.update()` | survives; types unchanged `['Ar0s', 'Ar+']` |
| `get_atomtype` per atom | `Ar0s` (neutral centre), `Ar+` (cationic centre) |
| `get_formula()` | `Ar2` |
| `get_net_charge()` | `+1` |
| `get_radical_count()` / `multiplicity` | `1` / `2` |
| `is_linear()` | `True` |
| `get_molecular_weight()` | `0.0797550 kg/mol` (= 2 × 39.878 g/mol) |
| `to_adjacency_list()` → `from_adjacency_list()` | round-trips; `is_isomorphic` **True** |
| `is_balanced()` on `Ar+ + Ar <=> Ar2+` | **True** |
| negative control `Ar+ <=> Ar2+ + e-` | **False** (so the True above is not vacuous) |

Round-tripped adjacency list (atom order normalised, isomorphic to the input):

```
multiplicity 2
1 Ar u1 p3 c0 {2,S}
2 Ar u0 p3 c+1 {1,S}
```

## 4. What still does not work

### 4a. SMILES corrupts the charge — the dimer **is** affected

```
to_smiles()  = '[Ar][Ar+]'
from_smiles  -> [('Ar', +1, 0, 3), ('Ar', +1, 0, 3)]   net charge +2   isomorphic back: False
```

Both argon atoms come back at `+1`, so the dimer cation reads back as a **dication**, and RDKit emits
`Explicit valence for atom # 0 Ar, 1, is greater than permitted` on the way. The monatomic control
run in the same script reproduces the known defect (`[Ar+]` → net `+2`), so the dimer is the *same*
defect, not a new one. **Not touched** — separate, still-gated ticket.

InChI, by contrast, is fine: `to_inchi()` gives `InChI=1S/Ar2/c1-2/q+1`, charge preserved.

### 4b. Thermochemistry — and the group-additivity question, answered

**Group additivity does NOT return a silent number for Ar₂⁺.** It refuses, three independent ways,
each loud:

1. `db.thermo.get_thermo_data_from_groups(Ar2+)` → `AtomTypeError`. It routes through HBI
   radical saturation, which caps the neutral radical centre with an H, producing a **two-bond**
   neutral argon (`Ar u0 p3 c0 {2,S} {3,S}`) — which `Ar0s` does not cover (`single=[0,1]`) and
   nothing else does either.
2. `db.thermo.get_thermo_data(Ar2+)` — the full library→QM→GA path RMG actually calls — fails at
   the same HBI step, same exception.
3. `db.thermo.compute_group_additivity_thermo(Ar2+)`, called directly to bypass HBI →
   `AssertionError: This method is only for saturated non-radical species.`

The deeper gate, established by the negative control on **monatomic** argon (which types fine and has
always typed fine, so the new atom type is not what stops it):

```
rmgpy.exceptions.DatabaseError: Unable to determine thermo parameters for atom {'*': <Atom 'Ar'>}
in molecule <Molecule "[Ar]">: no data for node R or any of its ancestors in database group.
```

**The thermo group tree contains no argon node at all.** So even if every atom-typing obstacle above
were removed, group additivity would still refuse rather than fabricate. The positive control in the
same run (`ethane`) returns `H298 = -20.400 kcal/mol` from
`group(Cs-CsHHH) + group(Cs-CsHHH)`, proving the estimator was live and the refusals are real.

This closes the ticket's loudest worry: **no**, nothing now silently invents a number for Ar₂⁺.

### 4c. …and no source publishes one

Stated plainly, as the brief asks: a buildable species is necessary but **not sufficient** for the
dimer electron-loss channel. No reachable source publishes a formation enthalpy for Ar₂⁺, and this
probe added none — no thermo entry, no database entry of any kind. Ar₂⁺ remains unusable in a
mechanism until thermochemistry for it exists, which is a data problem this change does not touch.

## 5. The test suites — and the one real finding in them

`pytest test/rmgpy/molecule/atomtypeTest.py test/rmgpy/molecule/moleculeTest.py -p no:cacheprovider --no-cov -q`

| | result |
|---|---|
| **before** (base `546d1c727`, rebuilt) | **233 passed, 3 skipped**, exit 0 |
| **after, atom type only** | **232 passed, 1 failed, 3 skipped**, exit 1 |
| **after, with the allowlist entry below** | **233 passed, 3 skipped**, exit 0 |

The one failure is worth more than the green that follows it.

`TestAtomType::test_make_sample_molecule` builds a sample molecule for every atom type from the group
`1 <name> ux` and asserts none crashes. For `Ar0s` it raises `UnexpectedChargeError`. The cause,
measured rather than guessed:

```
Ar0    OK  ->  1 Ar u0 p4 c0
Ar+    OK  ->  1 Ar u0 p3 c+1 {2,S} / 2 H u0 p0 c0 {1,S}
Ar++   OK  ->  1 Ar u0 p3 c+2
Ar0s   UnexpectedChargeError on:  1 Ar u0 p3 c+1 {2,S} / 2 H u0 p0 c0 {1,S}
```

`make_sample_molecule` reads the group's `ux` as `u0`, builds `Ar` with one bond and `p3`, and
`update_charge` then computes that atom to be **+1** — which the generator correctly refuses against
`Ar0s`'s declared `charge=[0]`.

**This is the guard reporting a true fact, not a false positive, and not a defect in the new type:
`Ar0s` has no closed-shell instance.** By the §2 arithmetic, a neutral argon holding one bond must
carry an unpaired electron; with zero bonds and `p3` it would be `+2`. There is no `u0` argon that is
both neutral and matches `Ar0s`. The generator's `ux → u0` assumption is what cannot be satisfied —
the same class of limitation as the `O4b`/`S4b` entries the test already carries. It was therefore
recorded in that same mechanism, with the reasoning in a comment at the site:

```python
EXPECTED_FAILING_ATOMTYPES = ["O4b", "S4b", "Ar0s"]
```

The forward-looking consequence, which is the part worth keeping: **anything in RMG that generates
sample molecules from atom types will hit this wall for `Ar0s`.** It is inert today only because no
database group names `Ar0s`: a `grep -rnE "^\s*[0-9]+\s+Ar(0|0s|\+|\+\+)\s"` over
`RMG-database-plasma/input` returns nothing — no adjacency-list line anywhere in the database uses an
argon atom-type label (the only byte match for `Ar0`/`Ar++` in the tree is inside the binary
`thermo/ml/main/s298_cp/model_3/model.pt`). It would stop being inert the moment a kinetics family
group did.

The atom-type table's own structural tests — `TestActionGraphClosure` (both-ways closure over the
whole table), `test_argon_and_alkali_families_close_both_ways`, `test_each_atomtype_label_defined_once`
— all pass unchanged, which is what the inert `set_actions` of §2 buys.

## 6. What the change widens — a full census, and it is two cells

Every argon local environment the table can be asked about was enumerated: bonds 0–2 × lone pairs 0–4
× charge 0..+2 × radicals 0–2 = **135** combinations, 14 of them buildable at all. Each was resolved
twice on the *same* loaded binary, differing in exactly one thing — `Ar0s` removed from, then present
in, `ATOMTYPES['Ar'].specific` (`widening-stdout.log`).

**Exactly two environments newly resolve, and both are electron-correct neutral argon:**

| bonds | p | c | u | before | after | own e⁻ | physical |
|---|---|---|---|---|---|---|---|
| 1 | 3 | 0 | 1 | *(unresolvable)* | **`Ar0s`** | 8 | ✔ — the Ar₂⁺ neutral half, the target |
| 0 | 3 | 0 | 2 | *(unresolvable)* | **`Ar0s`** | 8 | ✔ — a hypothetical neutral argon diradical |

**Zero regressions**: no environment that resolved before resolves to anything different now. `Ar0`,
`Ar+` and `Ar++` are untouched in every cell.

The second row is the whole of the collateral widening: a bond-free neutral argon carrying two
unpaired electrons would now type as `Ar0s` instead of erroring. It is electron-correct bookkeeping
and not a real species; whether admitting it is acceptable is a judgement for whoever decides on
keeping this change, not something the probe settles.

## 7. Side findings (not fixed — outside this ticket)

1. **`ConsistencyChecker.check_partial_charge` masks its own error.** At
   `rmgpy/molecule/adjlist.py:109` the `InvalidAdjacencyListError` message is built by calling
   `get_atomtype(atom, atom.edges)` — so when the atom's type is unknown, that call raises
   `AtomTypeError` and *replaces* the accurate valency diagnostic. The census caught this: four
   non-physical argon environments (6, 7 and 9 electrons) swap from `AtomTypeError` to the correct
   `InvalidAdjacencyListError` purely because `Ar0s` lets the message be constructed. Same rejection
   either way — but the pre-`Ar0s` diagnostic was actively misleading, and this masking is general,
   not argon-specific. It did **not** affect the §1 reproduction: that traceback comes from
   `molecule.py:1587`, not from the consistency checker.
2. **`RMGDatabase.load()` only loads thermo when `surface=True`** (`rmgpy/data/rmg.py:112`:
   `if surface: self.load_thermo(...)`). Calling `load(..., surface=False)` leaves `db.thermo` as
   `None` and the next thermo call dies with `AttributeError: 'NoneType' object has no attribute …`.
   The probe works around it by calling `load_thermo` directly.

## 8. What this probe may and may not be claimed to show

**May be claimed:**

- Ar₂⁺ was never structurally unrepresentable in RMG. The gap was one missing atom type for the
  bonded *neutral* argon; the bonded *cation* already typed, measured with a positive control.
- One atom type — `single=[0,1]`, `lone_pairs=[3]`, `charge=[0]` — makes Ar₂⁺ build, type, update,
  round-trip through its adjacency list, and balance in a reaction, on a freshly rebuilt binary
  against `/home/alon/Code/RMG-database-plasma/input`.
- Group additivity does not return a number for Ar₂⁺. It refuses three ways, and the underlying
  reason (no argon node anywhere in the thermo group tree) is broader than this species.
- The two named unit suites are green before and after, at identical counts, with the single
  difference explained and recorded rather than papered over.
- The change widens argon atom-typing by exactly two electron-correct neutral environments and
  regresses none.

**May NOT be claimed:**

- **That Ar₂⁺ is usable.** It has no thermochemistry, and none exists to attach. Buildable is
  necessary, not sufficient, for the dimer electron-loss channel.
- **That the dimer channel now works, or that any mechanism generates it.** No reaction, library,
  family, or deck was added or run. This probe never put Ar₂⁺ into a reactor.
- **That a family recipe can reach `Ar0s`.** Its `set_actions` are deliberately empty (§2), so no
  recipe moves an atom into or out of it. Only hand-written adjacency lists reach it.
- **That `Ar0s` is the right *design*.** Its name, its admission of the 0-bond diradical case, and
  the two undeclared inverse edges are all open choices. This probe measured one implementation; it
  did not compare it against alternatives.
- **That the wider RMG suite is unaffected.** Only `atomtypeTest.py` and `moleculeTest.py` were run.
  The atom-type tree is shared by every species in the engine; a green pair of unit files is the
  minimum bar, not evidence of safety at generation time. In particular the database-marked tests,
  the functional tests, and the regression models were **not** run.
- **Anything about `Ar₂⁺` beyond one adjacency list.** No resonance structures, no isomer search,
  no Chemkin/Cantera export, no transport, and no attempt at the three-body association that would
  actually form it.

**What the probe could not reach:** any statement about RMG's *runtime* behaviour with Ar₂⁺ present —
generation, resonance, export, or a solved reactor — because none of that was exercised, and 4b/4c
say why exercising it is blocked upstream on thermochemistry rather than on atom typing.
