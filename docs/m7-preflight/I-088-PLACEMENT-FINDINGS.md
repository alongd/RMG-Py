# I-088 — Electron placement for the cation family (Wall C)

**Date:** 2026-08-23 · **Branch:** `i088-electron-placement` · **Base:** `i086-electron-sign@0f5934dd2`
**Worktree:** `/home/alon/Code/RMG-Py-i088-placement` · **Database:** `/home/alon/Code/RMG-database-plasma` @ `fb3c13c60`

Scope: the reverse-generated lithium cation reaction `Li+ + CH3 (+ e-) <=> CH3Li` being refused by
`rmgpy/electron_placement.resolve_electron_placement` on three separate grounds. Evidence under
`repro-evidence/`. Regression lock: `test/rmgpy/electronPlacementTest.py`.

**Outcome: placement resolves.** The reaction now reaches the reactor as
`[Li+] + [CH3] + e- <=> CH3Li` with the electron on the reactant side, exactly once, and the view
balances in both net charge and the `E` pseudo-element. It still cannot be *integrated* — the
reactor refuses it for an independent and correct reason, named at the end of this note.

---

## Correction to the ticket's premise: the driver does not reach placement

The ticket states the driver fails at `plasma.pyx:324 -> electron_placement.py`, with a rate-law
export wall expected *behind* it. Measured at base with the committed deck
(`python rmg.py docs/m7-preflight/input.py`), that is **not what happens**. The run dies strictly
*earlier*, in the Chemkin **edge** write, on the Marcus rate-law wall:

```
main.py:961 save_everything -> chemkin.pyx:2004
MechanismWriterError: Cannot write reaction [Lip](3)+[CH3](4)+e-(1)<=>CH3Li(2) to Chemkin:
  Marcus kinetics ... dGrxn is not a property of the rate law
```
(`repro-evidence/I-088-BASE-deck-stops-at-marcus-writer-wall.stderr.log`)

That happens before the main loop's first `simulate()`, so the reactor never sees the reaction.
The one earlier `initialize_model` call (`main.py:912`, the reaction-filter pass) is invoked with
`edge_reactions=[]`, so the cation reaction cannot reach the resolver there either. **The expected
rate-law wall is in front of placement in the driver's execution order, not behind it.**

Placement *is* reachable from the driver with a one-flag **deck-config** change and no code change:
`saveEdgeSpecies=False` skips the edge Chemkin write (`chemkin.pyx:2645`), the run proceeds into the
main loop, and the ticket's failure appears verbatim:

```
main.py:1069 simulate -> base.pyx:678 -> plasma.pyx:237 -> plasma.pyx:334
  -> electron_placement.py:159
ElectronPlacementError: Family 'Cation_R_Recombination' has no electron-placement declaration
  (reaction [Lip](3) + [CH3](4) <=> CH3Li(2), electrons=-1)
```
(`repro-evidence/I-088-BEFORE-placement-wall.stderr.log`)

That deck is committed as `docs/m7-preflight/input_placement.py` — identical to `input.py` except for
the flag. Both decks are kept: `input.py` reproduces the Marcus writer wall, `input_placement.py`
reproduces (and now clears) the placement wall behind it.

---

## The three refusals, one at a time

### Refusal 1 — no `FAMILY_ELECTRON_PLACEMENT` entry (`electron_placement.py`, step 3)

**What it checks.** The family label has a `(side, count)` entry in the registry. Absence is a named
failure, never a fallback to inferring the side from `Reaction.electrons`.

**What it was protecting.** A real and correct invariant: the net scalar cannot distinguish
`A + e- -> A-` (one incident electron) from an ionization shape whose net production says the same
thing about the count but nothing about incident order. The registry is the order source.

**Disposition — declared.** `Cation_R_Recombination` is family-forward
`Li+ + R. + e- -> R-Li`, `electrons = -1` (`RMG-database-plasma`,
`kinetics/families/Cation_R_Recombination/groups.py:23`), i.e. one electron on the reactant side.
It gets `('reactants', 1)`, the same shape entry `Plasma_Electron_Attachment` already has. The
registry stays a closed, hand-maintained list — the other five cation families and the six
`Surface_Proton_Electron_Reduction_*` families all declare `electrons = -1` and remain undeclared,
locked by `test_declaration_registry_is_explicit_and_closed`.

This was the only refusal that was purely an omission. On its own it accomplishes nothing: refusals
2 and 3 stop the same reaction one and three lines later.

### Refusal 2 — `is_forward=False` (and `is_forward is not True`) (step 6)

**What it checks.** That the reaction was generated in the family's forward direction. Two raises: an
explicit `is False` refusal, and a stricter `is not True` refusal added later (BLOCKER 1) to catch
`is_forward=None`.

**What it was protecting — stated.** From the code comment: "a reverse-direction or reversible
reaction would put the electron on the wrong side, or on both", and "placing the electron on the
family-declared forward (reactant) side of a reaction whose direction was never established would
silently manufacture a forward-direction view from ambiguous input."

**What it was actually protecting — measured. Nothing that `is_forward` can see.** The premise —
that `is_forward=False` implies reversed participant lists — is false for every reaction that can
reach this resolver. All four producers of `is_forward = False` in the tree were read:

| producer | orientation stored | `electrons` | reaches resolver? |
|---|---|---|---|
| `family.py:1753` `_create_reaction` | **family-forward in both directions** (the reverse branch swaps the reverse-matched lists back) | family-forward `-1`, unnegated | **yes** — this is the engine's generation path |
| `family.py:4241`, `4283` `get_training_set(get_reverse=True)` | genuinely reversed (`reactants=products`) | negated to `+1` | no — plain `Reaction`, no family attribution (refused at step 1), and `+1` refused at step 5 |
| `family.py:1256` depository reverse entries | depository entries, not model reactions | n/a | no |

So `is_forward` records *how a match was found*, not how the lists are ordered, and it is
uninformative about which side the electron belongs on. Refusing on it excluded precisely the
chemistry the plasma decks generate — the family is `reversible = True` with an auto-derived reverse
(`reverse = "Cation_Bond_Dissociation"`), and it is the reverse generation direction that produces
the cation at all, so this refusal could never have been satisfied by this family.

Worse, the check does not even do the job it claims. A reaction that *is* genuinely
reverse-oriented while its metadata still claims the family-forward `-1` — the I-086 defect class —
passes every check including `is_forward=True`, and the base resolver **accepts** it and returns an
`E`-unbalanced view. Measured directly against `git show HEAD:rmgpy/electron_placement.py`:

```
BASE resolver ACCEPTED the wrong-side reaction -> O2- + e- => O2
   E left = 2   E right = 0
```

**Disposition — removed, and replaced with a stronger check.** What the `is_forward` refusals were
proxying for is "is the electron really a reactant in this reaction's stored orientation?" That is
now *verified* rather than trusted, at the new step 10: the finished view must balance in the `E`
pseudo-element, counted with the writers' own `get_species_electron_count` (imported read-only;
`electron_balance.py` is untouched). With *n* electrons on one side, moving them across changes the
imbalance by 2*n*, so an `E`-balanced view is proof the electron went to the correct side —
independent of `is_forward`, of `electrons`, and of whether the generating family honoured its
orientation invariant. The wrong-side reaction above is now refused by name; `is_forward=False` and
`is_forward=None` resolve. Both directions are locked
(`test_wrong_side_placement_refused_by_E_balance`, `test_reverse_generated_reaction_resolves`,
`test_unknown_direction_resolves`).

### Refusal 3 — `reversible=True` (step 6)

**What it checks.** That the reaction is irreversible.

**What it was protecting.** From the code comment: "a reversible view would leave the reverse rate of
an electron-containing reaction implicitly defined." **That hazard is real and it is not dismissed
here.** For a reversible reaction the reactor sets `kr = kf / Keq(Tgas)`; with the electron an
explicit participant, `Keq(Tgas)` prices the electron's thermochemistry at the *gas* temperature
while the electron population is held at `Te`. Two incompatible thermal closures in one rate.

**But it is not a placement question.** A reversible reaction still has a definite reactant side in
its stored orientation; which side the electron goes on is not ambiguous. Reversibility governs the
reverse *rate*, and that policy is owned and enforced by the consumer:

- `plasma.pyx:536-543` — `_validate_reactions` refuses a reversible electron-containing reaction by
  name with `NonEquilibriumReverseRateError`, with exactly the `Keq(Tgas)` reasoning above.
- `plasma.pyx:597-605` — `generate_rate_coefficients` re-checks it defensively.
- `plasma.pyx:501-512` — the same policy for `uses_electron_temperature` kinetics (RULING 1).

Both run on the **view**, after placement (`plasma.pyx:236-237` resolve, `:244` validate), and the
view carries `reversible` through unchanged, so placement cannot launder it. The resolver's copy of
the policy added nothing but an earlier refusal under a misattributed message: the declaration is
`('reactants', 1)` and says nothing about reversibility, yet the error read "family
'Plasma_Electron_Attachment' declares irreversible attachment".

**Disposition — removed as mislocated, with the protection verified still live.**
`TestCationRecombinationPlacement.test_05_reactor_still_refuses_the_reversible_view_by_name` takes
the real resolved view to a real `PlasmaReactor` and asserts `NonEquilibriumReverseRateError`. If the
reactor ever stops refusing, that test goes red — the protection is now locked at its owner rather
than duplicated at a boundary that does not own it.

---

## Where the pipeline stops now

`python rmg.py docs/m7-preflight/input_placement.py` gets the cation reaction into the reactor with
its electron correctly placed, and stops at the **reverse-rate wall**:

```
main.py:1069 simulate -> base.pyx:678 -> plasma.pyx:244 -> plasma.pyx:538
NonEquilibriumReverseRateError: automatic thermodynamic reversal is undefined for
  electron-containing reaction [Lip](3) + [CH3](4) + e-(1) <=> CH3Li(2), kinetics Marcus:
  Keq(Tgas) would price the electron's thermochemistry at the gas temperature;
  mark the reaction irreversible or provide explicit reverse kinetics
```
(`repro-evidence/I-088-AFTER-reverse-rate-wall.stderr.log`)

That equation is the deliverable: the electron is on the reactant side, once, and the reaction
balances. The wall behind it is not an electron-representation problem.

**What does not yet exist, precisely.** For this chemistry to be integrated in a two-temperature
plasma reactor, one of the following must be built, and neither exists today:

1. **Explicit reverse kinetics** for the electron-containing cation families, so no thermodynamic
   reversal is needed. This is what the error message asks for, and it is a database/chemistry
   deliverable, not a code one.
2. **A two-temperature equilibrium closure** — a `Keq(Tgas, Te)` defined for reactions whose
   participants are not all thermalised to one temperature — plus a ruling that it may be used. The
   current code takes the opposite position deliberately (RULING 1, `plasma.pyx:501`) and refuses
   the mixed closure rather than approximating it.

Marking these reactions irreversible is the third option the error names, but that is a chemistry
decision about the cation families and an explicit non-goal here.

Independently, the original deck (`input.py`, `saveEdgeSpecies=True`) still stops at the **Marcus
Chemkin wall**, unchanged and unchangeable from this side: the write happens before the reactor
runs, and Chemkin has no rate expression that reproduces Marcus's run-time `dGrxn` dependence. Two
distinct walls now remain, both out of scope, and they are reached by different decks.

---

## Durable finding (not fixed here): `Reaction.is_balanced` never compares charge

`rmgpy/reaction.py:1544-1596`. `is_balanced` declares and accumulates `reactants_net_charge` and
`products_net_charge` over every participant, then adjusts them by `self.electrons` — and
**returns `True` without ever comparing them**. Only the per-element counts are compared. So the
function's charge bookkeeping is dead code, and `is_balanced` cannot detect a charge imbalance.

This is why `_create_reaction`'s `if not reaction.is_balanced(): return None` did not catch the
I-086 electron-sign inversion, and why that defect had to be caught downstream at the Chemkin
writer. Not fixed here: making `is_balanced` compare charges would change reaction acceptance
across every family in the tree — a behavioural change far outside this ticket's Verifier, and one
that needs its own reproduction and blast-radius measurement.

---

## Verifier status

1. ✔ Failure reproduced at base **with the driver** (`input_placement.py`), plus the correction
   above: the committed `input.py` deck stops earlier, at the Marcus writer wall.
2. ✔ Per-refusal disposition with code evidence — the three sections above.
3. ✔ Placement resolves, locked by `TestCationRecombinationPlacement`, confirmed **RED before**
   (`repro-evidence/I-088-RED-before-placement-test.log`: 5 failed / 1 passed, the one pass being
   the precondition test that asserts only the reaction's shape) and green after.
4. ✔ Driver re-run; new stopping point named: `NonEquilibriumReverseRateError` at `plasma.pyx:538`.
5. ✔ `pytest test/rmgpy/electronPlacementTest.py test/rmgpy/solver/plasmaTest.py
   test/rmgpy/rmg/inputTest.py test/rmgpy/data/kinetics/` → **310 passed, 3 skipped**, exit 0, with
   `database.directory = /home/alon/Code/RMG-database-plasma/input`.
6. ✔ `git diff` on `rmgpy/electron_balance.py` and `rmgpy/chemkin.pyx` is empty. The E-balance guard
   is imported read-only and not modified, relaxed, or special-cased.
