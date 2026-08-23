# I-086 — The electron's sign on a reverse cation reaction (Wall D cleared)

**Date:** 2026-08-23 · **Branch:** `i086-electron-sign` · **Base:** `i085-compose-walls@8857c3516`
**Worktree:** `/home/alon/Code/RMG-Py-i086-electron-sign` · **Database:** `/home/alon/Code/RMG-database-plasma/input`

Scope: the third of the M7-preflight walls — the Chemkin writer's element-`E` balance guard
refusing the reverse-generated lithium cation reaction `[Li+] + [CH3] <=> CH3Li`. Reproduction and
before/after evidence under `repro-evidence/`. Regression lock:
`test/rmgpy/data/kinetics/cationReverseElectronSignTest.py`.

## The defect — a sign inversion in reaction generation, not the guard

**The guard was right; the representation was wrong.** At the writer, the reaction object was:

```
[Lip](3) + [CH3](4) <=> CH3Li(2)
  is_forward=False  family='Cation_R_Recombination'  electrons=+1  reversible=True
  reactant charges=[+1, 0]   product charges=[0]
```

`Cation_R_Recombination` is reversible and not its own reverse; its forward template is
`Li+ + R. (+ e-) -> R-Li` with the family-forward declaration `electrons = -1` (one electron
consumed, on the reactant side). With only the neutral product `CH3Li` in the model, RMG reaches
the family in the **reverse generation direction**: it matches `CH3Li` against the product template
and reconstructs the reactant side `Li+ + CH3`.

Crucially, `KineticsFamily._create_reaction` stores that reaction back in **family-forward
molecular orientation** — `reactants=[Li+, CH3], products=[CH3Li]` — while flagging it
`is_forward=False` to record *how it was found*. The empirical probe confirms the reactant side is
the family reactants (`[Li+, CH3]`), i.e. family-forward orientation. In that orientation the
equation only balances as `Li+ + CH3 + e- -> CH3Li`, so the electron is a **reactant** and
`Reaction.electrons` must be `-1`.

The bug: `_create_reaction` negated `electrons` whenever `is_forward` was false —
`electrons=self.electrons if is_forward else -self.electrons` — treating the `is_forward` flag as if
it also reversed the molecular orientation. It does not: the reactant/product lists were swapped
back to family-forward in the *same* constructor. The negation produced `electrons=+1` for a
reaction whose electron is physically a reactant. `expand_electrons` then places a `+1` electron on
the **product** side, giving `E=-1` on the left (from `Li+`) and `E=+1` on the right (from the
electron), and `check_electron_balance` correctly refused it.

### Why the negation is wrong in general (all real callers)

Every real caller invokes `_create_reaction(reactant_structures, product_structures, forward)` where
`reactant_structures` is what was matched against the forward *or reverse* template and
`product_structures` is the recipe output. For `forward=False` that is
`(family-products, family-reactants)`; the constructor's `reactants=... if is_forward else products`
swaps the reactant side back to the family-reactant molecules. So **both** directions store
family-forward orientation, and `Reaction.electrons` (signed to the stored orientation) is the
family-forward declaration `-1` in both. `_create_reaction` never produces a genuinely
reverse-oriented reaction.

The genuinely reverse-oriented reaction — reactant side = the *original products*, electron count
therefore negated to `+1` — is built elsewhere: the training-set flip in
`get_training_set(get_reverse=True)` at `family.py:4228/4269`, which constructs
`Reaction(reactants=products, ..., electrons=-rxns[i].electrons)` directly. **That path was not
touched** and its lock (`TestReverseTrainingReactionNegatesElectronCount`) still passes.

## The fix

One line, at the representation-generation site — **not** the guard:

- `rmgpy/data/kinetics/family.py`, `_create_reaction`: `electrons=self.electrons` (was
  `self.electrons if is_forward else -self.electrons`), with a comment recording the invariant and
  the I-086 defect. `family.py` is pure Python (not in `setup.py` `ext_modules`) — no rebuild.

`git diff` on `rmgpy/electron_balance.py` and `rmgpy/chemkin.pyx` is **empty**. The E-balance guard
now passes because the reaction balances, not because it changed.

### A correction to an existing test (it encoded the bug)

`TestCreateReactionNegatesElectronsWhenReversed` (in `electronPropagationTest.py`) called
`_create_reaction` with family-forward-oriented arguments *and* `is_forward=False`, then asserted
the negated `+1`. **No real caller uses that convention** — it inverts the `(matched-input,
recipe-output)` argument order the engine always uses, producing an artificial reverse-oriented
reaction. The test was renamed to `TestCreateReactionStoresFamilyForwardOrientation` and rewritten
to call `_create_reaction` the way the engine does and assert the true invariant (stored
family-forward, `electrons=-1`) in both directions. The correction is documented in the test's
docstring.

## The two quantities, asserted separately (the ticket's requirement)

`test/rmgpy/data/kinetics/cationReverseElectronSignTest.py` generates the reverse cation reaction
from the real family and asserts, by two independent routes:

- **Net electron change** — read off `Reaction.electrons`; must be `-1`.
- **Reactant-side electron participation** — obtained by expanding the electron into the writer's
  equation (`expand_electrons`) and counting: exactly `1` on the reactant side, `0` on the product
  side; corroborated by `get_molecularity == 3`.

Confirmed **RED before** the fix (all four quantity assertions failed: `electrons=+1`, electron on
product side, molecularity `2`, balance raised) and **GREEN after**.

## Wall C is a DIFFERENT defect (ticket hypothesis refuted)

The ticket hypothesised Wall C (`electron_placement.py`, reactor side) was "probably one missing
declaration observed at two sites." **It is not.** With the sign corrected (`electrons=-1`),
`resolve_electron_placement` on the same reaction still raises — at the *declaration* check:

```
ElectronPlacementError: Family 'Cation_R_Recombination' has no electron-placement declaration
  (reaction [Li+] + [CH3] <=> [Li][CH3], electrons=-1); ... Only families declared in
  FAMILY_ELECTRON_PLACEMENT can resolve.
```

Wall D was a **sign defect in reaction generation**; Wall C is a **missing family entry** in
`FAMILY_ELECTRON_PLACEMENT` (only `Plasma_Electron_Attachment` is declared). Moreover, even with a
declaration, `resolve_electron_placement` refuses this reaction on two further grounds by design:
it is `is_forward=False` and `reversible=True`. So clearing Wall C is a larger, separate piece of
work (a `Cation_R_Recombination` placement declaration plus reverse/reversible directionality
handling), not this ticket's one-line sign fix. The two walls are independent.

## Where the pipeline stops now — the Marcus wall (a new, out-of-scope wall)

After the fix the driver runs further: it generates the cation, the writer places the electron on
the reactant side, and the E-balance guard passes with **exact** balance —

```
[Lip](3)+[CH3](4)+e-(1)<=>CH3Li(2)      (E: 0 = 0, charge: 0 = 0)
```

— then stops at a **new** wall, unrelated to electron representation:

```
MechanismWriterError: Cannot write reaction [Lip]+[CH3]+e-<=>CH3Li to Chemkin: Marcus kinetics ...
  dGrxn is not a property of the rate law -- it comes from the species thermochemistry at run time.
  There is no value of these parameters for which a Chemkin rate expression reproduces this.
```

All `Cation_R_Recombination` training reactions carry **Marcus** kinetics (electron-transfer rate
law); Chemkin has no rate expression that reproduces Marcus's run-time `dGrxn` dependence, so the
writer correctly refuses it (`chemkin.pyx:2003`). This is a **rate-law-expressibility** limit, not
an electron-representation one, and the family's kinetics are an explicit non-goal of this ticket.
Before the fix, the E-balance wall masked it — the same "one wall masks the next" pattern the ticket
described for C and D.

## General "second free electron" limitation — NOT addressed

This fix corrects the **sign** of the net/reactant electron for the reverse recombination
(consumption) case. It does **not** give the reaction-family machinery the ability to emit a second
free electron (the ionization shape `A + e- -> A+ + 2 e-`), which `electron_placement.py` still
refuses as ionization-shaped and which the family machinery still cannot emit. An honest statement:
this fixes the reverse-family consumption case; the general second-free-electron limitation remains.

## Verifier status

1. ✔ Failure reproduced verbatim at base (`repro-evidence/I-086-BEFORE-fix-Ebalance-wall.log`).
2. ✔ RED-before regression test on the two separate quantities (RED confirmed, now GREEN).
3. ◑ Driver generates the cation; E-balance guard passes with exact charge/element balance and the
   electron on the reactant side. A written Chemkin mechanism is **not** reached — blocked by the
   independent, out-of-scope Marcus wall, not by electron representation.
4. ✔ `git diff` on `electron_balance.py` and `chemkin.pyx` is empty; the fix is one line in
   `family.py`.
5. ✔ `pytest test/rmgpy/solver/plasmaTest.py test/rmgpy/rmg/inputTest.py` green (with
   `database.directory` → RMG-database-plasma); full `test/rmgpy/data/kinetics/` +
   `reactionTest.py` also green (215 passed, 4 skipped).
6. ✔ Wall C: does NOT clear with this fix; it is a separate missing-declaration defect (shown
   above). General second-free-electron limitation: remains.
7. N/A — the export does not succeed (Marcus wall), so the downstream load with the repair
   workaround disabled cannot be reached from this deck. The electron representation no longer
   blocks export; the Marcus rate law does, independently.
