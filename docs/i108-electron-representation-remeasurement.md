# I-108 — Electron placement and ionisation: re-measurement

Base: `i108-electron-placement` @ `0a3b0ff3d` ("Let the two plasma rate laws declare the electron
they consume"). Worktree `/home/alon/Code/RMG-Py-i108-electron-placement`, database
`/home/alon/Code/RMG-database-plasma/input`, `rmg_env`.

Two investigations measured the same machinery and drew opposite conclusions; one measured against
a base that has since moved. This report re-measures rather than adjudicates. Every claim is
labelled **[R]** code, **[D]** database, **[M]** measured, **[I]** inference.

The measurement harness is committed at
[`test/rmgpy/i108ElectronRepresentationMatrixTest.py`](../test/rmgpy/i108ElectronRepresentationMatrixTest.py).
It runs standalone (prints the tables below) and under pytest (asserts them).

---

## 0. Environment and build

**[M]**

```
$ python utilities.py check-pydas && make build
BUILD_EXIT=0

$ python -c "import rmgpy; print(rmgpy.__file__)"
/home/alon/Code/RMG-Py-i108-electron-placement/rmgpy/__init__.py

$ python -c "from rmgpy import settings; print(settings['database.directory'])"
/home/alon/Code/RMG-database-plasma/input
```

The worktree had **no** `.so` files before this (`find . -name '*.so' | wc -l` → `0`); after
`make build`, 104. Every number below was taken against that build. `rmgpy/reaction.py`,
`rmgpy/kinetics/arrhenius.pyx` and `rmgpy/solver/plasma.pyx` are all compiled, so this matters.

---

## 1. The green suite — corroborated, with a +1 that is accounted for

**[M]** On the base, `pytest test/rmgpy/` (default `pytest.ini` flags, no `-p no:cov`, my own added
test file not yet present at collection time — verified: `grep -c i108ElectronRepresentationMatrix`
on the log → `0`):

```
collected 2510 items
========== 2477 passed, 33 skipped, 57 warnings in 354.35s (0:05:54) ===========
```

The shipped change claims `2476 passed, 0 failed, 0 errors`. **We measured 2477.** Per the ticket
that discrepancy outranks everything, so it was chased before anything else was done.

**[M]** I rebuilt and re-ran the suite at the parent commit `6cb64c70d` — the base the shipped
change measured its "before" against:

```
$ git switch --detach 6cb64c70d && make build && pytest test/rmgpy/
===== 2469 passed, 34 skipped, 57 warnings, 5 errors in 344.11s (0:05:44) ======
```

That is **exactly** the claimed before-figure, `2469 passed + 5 errors`. The "before" reproduces
to the test.

**[M]** Diffing the two skip lists gives a single entry, present at the parent and absent on the
base:

```
$ comm -23 skips-parent.txt skips-base.txt
test/rmgpy/data/thermoTest.py::TestMolecularManipulationInvolvedInThermoEstimation::test_deterministic_bicyclic_decomposition
```

**[R]** That test self-skips nondeterministically. `test/rmgpy/data/thermoTest.py:2453-2477`:

```python
        # Ensure that the order is the same every time
        try:
            assert ring_smiles == ["C1C=CC=C=C1", "C1C=CCC=C1"]
        except AssertionError as e:
            pytest.skip(f"Skipping because not yet deterministic (#2562): {e}")
```

Its own docstring says *"Currently this is not guaranteed, so if this test fails, we just skip it"*
(RMG-Py issue #2562).

**Verdict on Task 1 — the claim reproduces, and the +1 is that flaky test.** `2469 + 5 rescued +
2 new = 2476`; we saw 2477 because this test happened to land on the deterministic ring order in
that run and pass rather than skip. It touches ring decomposition in thermo estimation, not
electrons.

**[M]** Confirmed directly by re-running the suite a second time on the same build, with this
ticket's 22 added tests present:

```
========== 2498 passed, 34 skipped, 57 warnings in 349.81s (0:05:49) ===========
```

`2498 - 22 = 2476` passed with `34` skipped — **exactly the shipped figure**, on a run where the
flaky test skipped instead of passing. Two runs of the same build, same commit, differing by one
test that flips between pass and skip; the claimed number is one of the two values it produces.
Failures: 0. Errors: 0, in both runs. The green suite is corroborated.

---

## 2. The matrix, re-measured on this base

`python test/rmgpy/i108ElectronRepresentationMatrixTest.py`. The three questions are asked of three
genuinely independent gates: `Reaction.is_balanced()`; `KineticsLibrary.load` (writes a one-entry
library to a temp dir and loads it); `PlasmaReactor.initialize_model`.

**[D]** The ionisation and recombination rows use real shipped lithium kinetics, not placeholders:
`VoronovEIArrhenius(Z=3, N=3)` reads the Li I → Li II fit from
`RMG-database-plasma/input/kinetics/voronov.yaml` (`dE 5.4 eV, A 0.139e-6, X 0.438, K 0.41`), and
`BadnellRRArrhenius(Z=3, N=2)` reads the Li II → Li I fit from `input/kinetics/badnell.yaml`
(`A 9.349e-11, B 0.6916, T0 34.70, T1 7.329e6, C 3.51e-2, T2 6.293e5`). Both resolve their path
through `settings['database.directory']` **[R]** (`arrhenius.pyx:887`, `arrhenius.pyx:1541`).

### 2.1 The eight rows **[M]**

```
     form             equation                                 is_balanced   library   reactor
----------------------------------------------------------------------------------------------
EI   explicit         Li + e => Li+ + e + e                         refuse    refuse    ACCEPT
EI   half             Li => Li+ + e                                 refuse    refuse    refuse
EI   metadata         Li => Li+  (electrons=+1)                     ACCEPT    ACCEPT    refuse
RR   explicit         Li+ + e => Li                                 refuse    refuse    ACCEPT
RR   metadata         Li+ => Li  (electrons=-1)                     ACCEPT    ACCEPT    refuse
3b   e-e-ion          Li+ + e + e => Li + e                         refuse    refuse    ACCEPT
3b   neutral M        Li+ + e + Ar => Li + Ar                       refuse    refuse    ACCEPT
CT   ion-atom         Li+ + Ar => Li + Ar                           refuse    refuse    ACCEPT
```

Per-row detail, with the reason for each refusal:

| row | `is_balanced` | why | library | reactor |
|---|---|---|---|---|
| EI explicit | refuse | element mismatch on `e` (L 1 / R 2) | refuse — *"not balanced! Please reformulate."* | **ACCEPT** |
| EI half | refuse | element mismatch on `e` (L 0 / R 1) | refuse — same | refuse — `PlasmaStateError`: Voronov *"declares `uses_electron_density`, but no explicit electron appears among its reactants"* |
| EI metadata | ACCEPT | balanced | **ACCEPT** — stored `electrons=+1` | refuse — `ElectronPlacementError`: *"carries no family attribution"* |
| RR explicit | refuse | element mismatch on `e` (L 1 / R 0) | refuse — same | **ACCEPT** |
| RR metadata | ACCEPT | balanced | **ACCEPT** — stored `electrons=-1` | refuse — `ElectronPlacementError`: *"carries no family attribution"* |
| 3b e-e-ion | refuse | element mismatch on `e` (L 2 / R 1) | refuse — same | **ACCEPT** |
| 3b neutral M | refuse | element mismatch on `e` (L 1 / R 0) | refuse — same | **ACCEPT** |
| CT ion-atom | refuse | **charge** mismatch (L +1 / R +0) — no element mismatch at all | refuse — same | **ACCEPT** |

Two things in that table were not in the earlier report and are worth stating plainly.

**[M] The two three-body rows are not electron-conserving.** `Li+ + e + e => Li + e` goes 2 free
electrons → 1, and `Li+ + e + Ar => Li + Ar` goes 1 → 0. They are refused by `is_balanced` for
exactly the same reason as `RR explicit`, not for a different one. There is no shape in the matrix
in which an explicit electron merely passes through.

**[M] The `CT ion-atom` row as written does not conserve charge.** `Li+ + Ar => Li + Ar` is +1 on
the left and 0 on the right, and contains no electron anywhere. Its refusal is a charge-balance
result and says nothing about electron representation. The charge-conserving spelling of that same
chemistry, `Li+ + Ar => Li + Ar+`, is **accepted by all three gates** (supplementary table below).
Reading row 8's `refuse` as evidence about the electron machinery would be a misreading.

### 2.2 Supplementary rows **[M]**

```
     form             equation                                 is_balanced   library   reactor
----------------------------------------------------------------------------------------------
CT   charge-balanced  Li+ + Ar => Li + Ar+                          ACCEPT    ACCEPT    ACCEPT
EI   half-explicit    Li + e => Li+ + e  (electrons=+1)             ACCEPT    ACCEPT    refuse
3b   e-e-ion+meta     Li+ + e + e => Li + e  (electrons=-1)         refuse    refuse    refuse
```

`3b/e-e-ion+meta` is the demonstration that **the two representations exclude each other rather
than compose**. **[R]** `Reaction.is_balanced` (`rmgpy/reaction.py:1622-1638`) compares per-element
atom counts *first* and returns `False` immediately; the signed scalar is folded only into the
*charge* comparison afterwards. So a shape that already fails on element `e` is beyond the scalar's
reach — measured detail: `element mismatch on e (L 2 / R 1); charge mismatch (L -2 / R -1 after
folding electrons=-1)`. Adding the scalar made it worse, not better.

`EI/half-explicit` is the sharpest row in the whole ticket and is discussed in §5.

### 2.3 What changed relative to the earlier measurement, and what did not

**Unchanged [M]:**

- `is_balanced()` still refuses every explicit-electron form whose free-electron count changes, and
  still refuses it on the **element** comparison, before charge is looked at.
- `is_balanced()` still accepts both metadata forms (`Li => Li+` with `+1`, `Li+ => Li` with `-1`).
- The reactor still refuses every nonzero metadata count, by name.
- `FAMILY_ELECTRON_PLACEMENT` still holds exactly two entries, both `('reactants', 1)`.

**Changed [M]:** the **library column for the two metadata rows**. The earlier report found *"no
metadata escape on the library route"*. On this base there is one, and both metadata rows load.
This is the one cell the shipped change moved. §3 settles it.

**Never measured before, and it changes the shape of the problem [M]:** the reactor column
**ACCEPTs** five of the eight rows — including `EI explicit` and `RR explicit`, the two forms
`is_balanced()` refuses. `PlasmaReactor` never calls `is_balanced()`; grep confirms it **[R]**. The
two gates are orthogonal in both directions: the reactor accepts what the balance check refuses,
and refuses (`RR metadata`) what the balance check accepts. Any account of this machinery as a
single "pincer" with two jaws that both have to open on the same representation is wrong on the
measurement — the representation each jaw wants is *different*, and neither is a superset of the
other.

---

## 3. Is the library route still shut? **No.** **[M]**

The earlier finding was: *"`KineticsLibrary.load_entry` has no `electrons` keyword, so a library
entry can carry the electron only as an explicit species, the one form `is_balanced()` refuses."*

The first half is still literally true **[R]** — `load_entry` (`rmgpy/data/kinetics/library.py:602`)
takes no `electrons` argument, and there is no `electrons=` keyword on the `entry()` DSL. The
conclusion no longer follows, because the count arrives by a different channel.

**[R]** `KineticsLibrary.load`, `rmgpy/data/kinetics/library.py:567-583`, immediately *before* the
balance check:

```python
            if hasattr(entry.data, 'electrons'):
                ...
                rxn.electrons = entry.data.electrons.value

            if not rxn.is_balanced():
                raise DatabaseError('Reaction {0} in kinetics library {1} was not balanced! ...')
```

The escape is borne by the **rate law**, not by the entry signature. Commit `04bbb2881`
("Propagate a library entry's declared electron count to its reaction") opened that channel and
commit `0a3b0ff3d` gave `BadnellRRArrhenius` and `VoronovEIArrhenius` the `electrons` field that
feeds it.

**[M]** Measured by writing a one-entry library carrying the real Li kinetics into a temp dir and
loading it:

```
== EI/metadata: Li => Li+  (electrons=+1)
   library     : ACCEPT  stored electrons=+1 on Reaction    # kinetics = VoronovEIArrhenius(Z=3, N=3)

== RR/metadata: Li+ => Li  (electrons=-1)
   library     : ACCEPT  stored electrons=-1 on Reaction    # kinetics = BadnellRRArrhenius(Z=3, N=2)
```

**[M]** And that it really is the rate law doing the work — the same equation with a plain
`Arrhenius`, which has no `electrons` attribute, loses the count and is refused by the loader's own
`is_balanced()` call:

```
Reaction Li => Liplus in kinetics library i108-probe was not balanced! Please reformulate.
```

(`test_library_metadata_escape_is_rate_law_borne`.) Note also that the entry is reconstructed from
`repr(entry.data)` by the loader's `exec`, so this round-trips through the persisted form — which
is why `0a3b0ff3d` made `repr` emit the count unconditionally.

**Answer to Verifier 4: yes.** The kinetics-library route can store an electron-count-changing
reaction on this base, demonstrated by loading two of them built on shipped database kinetics.
**The loader jaw is open.** **[I]** The limitation that remains is scope, not principle: only a
rate law that carries an `electrons` attribute can deliver the count. **[M]** Enumerated over
`rmgpy.kinetics`, exactly six classes qualify today:

```
['ArrheniusChargeTransfer', 'ArrheniusChargeTransferBM', 'BadnellRRArrhenius',
 'SurfaceChargeTransfer', 'SurfaceChargeTransferBEP', 'VoronovEIArrhenius']
```

and **[R]** wrapper kinetics such as `MultiArrhenius` expose no top-level `electrons` and are
deliberately not traversed (comment at `library.py:576-579`). **[R]** The same escape exists on the
depository route (`rmgpy/data/kinetics/depository.py:224-234`), which additionally falls back on
the owning family's declared count; a library has no owning family to ask.

---

## 4. The reactor jaw — untouched, all three refusals confirmed **[M]**

**[M]** `FAMILY_ELECTRON_PLACEMENT` is byte-for-byte what the ticket describes, asserted by
`test_family_electron_placement_is_still_the_exact_declared_table`:

```
FAMILY_ELECTRON_PLACEMENT = {'Plasma_Electron_Attachment': ('reactants', 1),
                             'Cation_R_Recombination': ('reactants', 1)}
```

The value above is **as of I-108** and has since moved twice, both times by hand: I-113 respelled it
to the two-sided `(reactant_count, product_count)` form, and I-116 added a third entry,
`'PlasmaElectronImpactIonization': (1, 2)` — a kinetics library, and the first electron-producing
owner. The test carrying this row was renamed at I-116 because its old name asserted "exactly two
consuming entries", which had become false. What the row still measures — the table is closed and
hand-maintained — is unchanged.

**[M]** All three refusals fire, with the messages intact (full output in the harness's `THE
REFUSALS` section):

1. **Undeclared family** — `Li => Li+`, `electrons=+1`, family `Li_Electron_Impact_Ionization`:
   > `Family 'Li_Electron_Impact_Ionization' has no electron-placement declaration (reaction Li => Liplus, electrons=1); refusing to infer electron placement from the net electron count.`

2. **Ionisation-shaped against a consumption-declaring family** — same reaction, family
   `Plasma_Electron_Attachment`:
   > `Reaction Li => Liplus carries electrons=1, which is ionization-shaped (net electron production); family 'Plasma_Electron_Attachment' declares single-electron consumption (net -1). No placement view is defined for this shape in this increment.`

3. **Reactor metadata-only count** — `Li+ => Li`, `electrons=-1`:
   > `PlasmaStateError: reaction Liplus => Li carries a metadata-only electron count (electrons=-1); this representation cannot distinguish incident-electron order from net electron production, so it is unsupported here.`

**[R] One qualification on refusal 3.** It is reached in the harness by calling
`PlasmaReactor._validate_reactions` directly. On the whole-pipeline path it is *unreachable*:
`initialize_model` calls `_resolve_electron_placements` first (`plasma.pyx:236-237`), which routes
every `electrons != 0` reaction to `resolve_electron_placement`, and that either raises or returns
a view with `electrons = 0`. So refusal 3 is a backstop behind the resolver, not the first thing a
metadata reaction meets — which is why the matrix's `EI metadata` and `RR metadata` rows report the
resolver's *"carries no family attribution"* message and not this one. That is correct defence in
depth, but it means a report that quotes refusal 3 as "what the reactor says to a metadata
reaction" is quoting the second line of defence.

**[M] A fourth refusal, not in the ticket's list, found by measuring.** See §5.

---

## 5. The finding: the export boundary prescribes what the placement boundary refuses

**[R]** `rmgpy/electron_balance.py:203-216`, `check_electron_reactant_order`, states the correct
RMG representation of electron-impact ionisation in its own docstring:

> *"The correct RMG representation puts the consumed electron in `reaction.reactants` and counts
> only the surplus produced electrons in `reaction.electrons` (`A + e => A+ + 2 e` is
> `reactants=[A, e]` with `electrons=2`)."*

**[D]** The shipped fixture agrees — `test/.../plasma-local-context/reactions.py` entry 4 is
`e + H => proton + e` with `electrons = 1`, and `0a3b0ff3d`'s message calls writing that incident
electron out *"load-bearing rather than decorative — it is what makes the equation second order,
as the cm^3/(molecule*s) rate coefficient requires."*

**[M]** That form balances and loads:

```
EI   half-explicit    Li + e => Li+ + e  (electrons=+1)             ACCEPT    ACCEPT    refuse
```

**[M]** And the placement resolver refuses it outright — with a declared family, so that step 1's
family check cannot mask what step 2 does:

```
Reaction Li + e => Liplus + e (family 'Plasma_Electron_Attachment') represents the electron twice:
2 explicit electron participant(s) AND a nonzero metadata electron count (electrons=1). This
double representation would double-count electron stoichiometry or rate order; refusing to prefer
either source.
```

**[R]** `rmgpy/electron_placement.py` step 2. And **[R]** `PlasmaReactor._validate_reactions`
refuses any `rxn.electrons != 0` besides.

**[I]** So the one representation that carries *both* pieces of information the physics needs —
incident order as an explicit reactant, net production as the scalar — is loadable, exportable,
and endorsed by two docstrings and a shipped fixture, and is **structurally unreachable through the
plasma reactor**, at any family, today. It is not that no family declares it; it is refused before
any declaration is consulted. This is a fourth jaw and, on the measurement, the most consequential
one, because it is the jaw that closes on the *right* answer.

### 5.1 Why the two prior investigations disagreed: there are two "E" elements **[M]**

The codebase carries two incompatible definitions of what balancing the electron means.

**[R]** `Reaction.is_balanced` counts **atoms of element `e`** (`rmgpy/reaction.py:1595`,
`1622-1624`); only the electron pseudo-species has any. Under this rule a free electron is a
conserved element: it may pass through but cannot be created or destroyed.

**[R]** `rmgpy.electron_balance.get_species_electron_count` (`electron_balance.py:113-128`) returns
**minus the net charge** for a charged species, falling back to counting `E` atoms only for neutral
ones. This is the rule used by `check_electron_balance` at the Chemkin/Cantera export boundary
**and** by `resolve_electron_placement` step 10 — the resolver's own structural verification that
the declared side was the right side.

**[M]** They disagree on every explicit-electron shape in the matrix
(`test_the_two_E_rules_disagree_on_explicit_ionisation`), `(left, right)` under each rule:

| shape | atom rule (`is_balanced`) | charge rule (writers, resolver step 10) |
|---|---|---|
| `Li + e => Li+ + e + e` | (1, 2) — **refuse** | (1, 1) — **balanced** |
| `Li => Li+ + e` | (0, 1) — **refuse** | (0, 0) — **balanced** |
| `Li+ + e => Li` | (1, 0) — **refuse** | (0, 0) — **balanced** |
| `Li+ + e + e => Li + e` | (2, 1) — **refuse** | (1, 1) — **balanced** |
| `Li+ + e + Ar => Li + Ar` | (1, 0) — **refuse** | (0, 0) — **balanced** |

**[I]** This is the whole of the disagreement. The first investigation looked at the explicit form
and saw a chemically correct, E-balanced reaction being refused, and called `is_balanced` the
defect. The second looked at `is_balanced`'s own rule and saw a conserved element being correctly
conserved, and called it load-bearing. **Both are right about the rule they were looking at.**
Neither is a defect in isolation; what is real is that the two boundaries do not agree on what the
electron *is*, and the placement resolver — the piece that exists precisely to translate between
representations — verifies its output with the charge rule while the database gate that produced
its input used the atom rule.

**Recommendation, not an edit** (the ticket forbids touching `is_balanced`, and this report does
not): the divergence should be documented as a deliberate two-rule design, with each rule's domain
named at both sites, or the two should be reconciled. It should not be left implicit, because every
future reader will re-derive one of the two prior conclusions from whichever site they read first.

---

## 6. Design recommendation — a placement view for ionisation-shaped reactions

**Recommendation only. Nothing below is implemented; `FAMILY_ELECTRON_PLACEMENT` and
`resolve_electron_placement` are untouched on this branch.**

### The problem, precisely

**[R]** The declaration is `(side, count)` and step 5 validates `reaction.electrons == -count`.
Attachment is the case where net change and incident order coincide: `A + e- -> A-` has one
electron on the reactant side, that is the whole of the net change, and it is also the whole of the
electron's contribution to reaction order. One number does both jobs, so the two-field shape has
never had to distinguish them.

Ionisation separates them. `Li + e- -> Li+ + 2e-` has **incident order 1** (which fixes the rate as
second order overall, and is what makes the `cm^3/(molecule*s)` Voronov coefficient dimensionally
meaningful **[D]**) and **net production +1**. A bare `electrons = +1` cannot tell that apart from
`Li -> Li+ + e-`, which is first order and would need units of `s^-1`. The rate would then be wrong
by a factor of the electron density while the file looks well formed — **[R]** exactly the failure
`check_electron_reactant_order` was written to catch at export.

### Alternatives

**A. Two-sided count: replace `(side, count)` with `(reactant_count, product_count)`. — recommended.**

- `Plasma_Electron_Attachment: (1, 0)`; `Cation_R_Recombination: (1, 0)`; an ionisation family:
  `(1, 2)`.
- Incident order is `reactant_count`, *declared* rather than derived. Net is
  `product_count - reactant_count`, and step 5's validation generalises from
  `electrons == -count` to `electrons == product_count - reactant_count` with no change of
  character: the count is still only ever *validated against* the declaration, never the order
  source. That property is the module's central invariant **[R]** and this preserves it exactly.
- **[M]** Step 10's structural verification generalises with **no edit at all**. Under the charge
  rule, `Li + e- -> Li+ + 2e-` is E=(1, 1) — balanced — so the finished ionisation view passes the
  same check the attachment view passes today. This was measured, not assumed
  (`e_counts_by_both_rules`), and it matters: step 10 is the check that would have caught the I-086
  sign inversion, and an extension that broke it would be trading the strongest guarantee in the
  module for a feature.
- **Cost.** The real one is **[R]** step 9's list-append primitive. `expand_electrons` is
  net-derived — it appends `-electrons` electrons to the reactants and nothing to the products —
  and it is shared with the export path, so it cannot express a two-sided placement. This needs
  either a second primitive or a widened one, and widening a function the writers also call is the
  risk to weigh. Beyond that: a schema change to a two-entry table (2 rows today), the shape check
  at step 3 becomes a general one, and every call site of the declaration.

**B. Add a third field: `(side, count, incident)`.**

- Smaller diff. But for ionisation you would write `('products', 2, 1)` — "two on the product side,
  one of which is also on the reactant side" — and the fields stop being independent: `count`,
  `incident` and `side` now have to satisfy an invariant that a reader must know and a validator
  must enforce. It encodes the same information as A with an extra consistency rule to maintain.
  **Rejected.**

**C. Derive incident order from the rate coefficient's units.**

- **[R]** `get_plasma_rate_order(kinetics)` already does exactly this, and
  `check_electron_reactant_order` already trusts it at export.
- As the *order source* this is wrong: it makes placement kinetics-derived, which the module's
  docstring explicitly forbids, and **[R]** it returns `None` for any rate law whose units are not
  in `_ORDER_BY_RATE_UNITS` — a silent gap, not a named failure.
- As a **cross-check** it is right, and should be adopted alongside A: after building the view,
  assert `get_plasma_rate_order(kinetics) == len(view.reactants)` where the order is resolvable,
  and refuse by name where it disagrees. That is the same shape the module already uses three
  times over — declare, validate the metadata against the declaration, then verify the finished
  object against the participants. It converts the export-time `MechanismWriterError` into a
  placement-time refusal, which is strictly earlier and cheaper.

**D. Accept the half-explicit canonical form (incident electron explicit, surplus as scalar) and
place only the surplus.**

- This is the form the export boundary prescribes and the fixture writes (§5), so it needs no
  database migration at all — the data is already written this way.
- **[R]** But it collides head-on with step 2, which refuses explicit-plus-scalar as double
  representation, and that refusal is protecting something real: with an arbitrary family, "2
  explicit electrons and `electrons=+1`" genuinely is ambiguous about whether the scalar is a
  surplus or a duplicate.
- **[I]** It becomes safe **only if** the declaration says which reading applies — i.e. only on
  top of A. With `(1, 2)` declared, "one explicit reactant electron plus `electrons=+1`" has
  exactly one consistent reading and step 2 can narrow from "any explicit electron with a nonzero
  count is fatal" to "any explicit electron the declaration does not account for is fatal".
- **Recommend as a follow-on to A, not an alternative to it.** It is the difference between the
  extension being usable on the data that exists and requiring every ionisation entry to be
  rewritten into the fully-collapsed metadata form.

### How A handles the incident-electron-order problem

The declaration carries it directly and separately: `reactant_count` is the incident order,
`product_count - reactant_count` is the net change, and neither is inferred from the other. The
scalar `Reaction.electrons` stays what it is — a net count, validated against the declaration,
never an order source. The finished view is then checked twice against things no producer
controls: the E pseudo-element balance (step 10, unchanged) and the rate coefficient's own
dimensionality (option C). **[I]** Three independent sources — family declaration, database scalar,
rate-law units — that must agree, which is the design the module already has, extended to the one
axis it currently cannot express.

**[I] The existing two-field shape cannot express it.** `(side, count)` names one side and one
number; ionisation needs a number on each side, and they are not derivable from one another. This
is a widening, not a reinterpretation.

### What this does not settle — UNKNOWN

- **Which family, and whether one exists.** No ionisation family is declared anywhere; the
  measurement used a synthetic label. **[D]** Whether RMG-database-plasma has, or should grow, an
  electron-impact-ionisation family — as opposed to library entries only — was not established
  here. What would settle it: an inventory of the plasma database's family set against the decks
  that need Li ionisation.
- **Whether a library reaction can ever resolve.** **[R]** `LibraryReaction` sets
  `self.family = library`, i.e. the library's label, and `resolve_electron_placement` looks that up
  in `FAMILY_ELECTRON_PLACEMENT`. So a library entry would resolve only if a *library label* were
  added to a family-keyed table. Whether that is intended, or whether library reactions are meant
  to reach the reactor already in explicit form, is not decided by anything in the code. What would
  settle it: the intended provenance of the plasma decks' ionisation reactions — family-generated
  or library-supplied.
- **Multi-electron ionisation.** `Li+ -> Li++` steps exist in `voronov.yaml` (N=2, N=1) **[D]**.
  A declares `(1, 2)` per step, which covers them one charge state at a time; nothing here checks
  that a family can be declared per charge state rather than per element.

---

## 7. Verdict

**Yes — a defensible alkali cation source is blocked by representation on this base, and the jaw is
the reactor, not the loader.** The loader jaw, which the earlier investigation found shut, is now
open: measured here, a kinetics library entry carrying real `voronov.yaml` electron-impact
ionisation stores `electrons=+1` and one carrying real `badnell.yaml` radiative recombination
stores `electrons=-1`, both surviving `is_balanced()` (§3). What blocks the cation source is
downstream of that, and it is not one refusal but a closed ring: the fully-collapsed metadata form
`Li => Li+ (electrons=+1)` loads but the reactor refuses every nonzero metadata count and no
ionisation family has a placement declaration to resolve it with; the fully-explicit form
`Li + e => Li+ + e + e` is what the reactor actually wants — measured, it *initializes* — but no
loader will produce it because `is_balanced()` refuses it on the element-`e` comparison; and the
half-explicit form `Li + e => Li+ + e (electrons=+1)`, which is the one the export boundary's own
docstring and the shipped fixture both call correct, loads perfectly and is then refused by the
placement resolver as double representation before any family declaration is even consulted. Three
representations, three different gates, and each representation is refused by a different one.
**[I]** The unblock is therefore not a single edit to any one gate: it is a placement declaration
that can express incident order separately from net change (§6, option A), which is what lets the
reactor be handed the explicit form it already accepts without asking the database to store a shape
`is_balanced()` will not pass. `is_balanced` needs no change for that, and this report recommends
none.
