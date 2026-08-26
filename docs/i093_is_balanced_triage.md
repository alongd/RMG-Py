# I-093 — `Reaction.is_balanced` charge check: fix + fallout triage

## The fix

`rmgpy/reaction.py::Reaction.is_balanced` computed a net-charge balance (with an electron-count
adjustment) and then discarded it with an unconditional `return True`. The charge half of the
function was dead: every reaction that passed the element check was reported balanced regardless
of charge. The fix is one line — `return reactants_net_charge == products_net_charge` — plus a
comment stating the electron-adjustment sign convention.

### Sign convention (stated and justified)

`Reaction.electrons` is the stoichiometric coefficient of the electron, signed to the object's
*current* orientation: **negative = electron consumed** (reactant side), **positive = electron
produced** (product side). This matches `rmgpy.electron_balance.expand_electrons`. An electron
carries charge −1, so |electrons| electrons on the reactant side contribute total charge
`self.electrons` (a negative number), and |electrons| on the product side contribute
`-self.electrons`. Hence the pre-existing adjustment:

```python
if self.electrons < 0:   reactants_net_charge += self.electrons     # consumed
elif self.electrons > 0: products_net_charge  -= self.electrons     # produced
return reactants_net_charge == products_net_charge
```

**One consuming reaction** — A + e⁻ → A⁻ (`electrons = -1`): reactant species charge 0,
`+= -1` → −1; product species charge −1; −1 == −1 → **balanced**.

**One producing reaction** — A⁻ → A + e⁻ (`electrons = +1`): reactant species charge −1;
product species charge 0, `-= +1` → −1; −1 == −1 → **balanced**.
(Equivalently A → A⁺ + e⁻: reactant 0; product +1, `-= 1` → 0; 0 == 0 → balanced.)

## RED-before evidence

`test/rmgpy/reactionTest.py::TestIsBalancedChargeConservation` — five tests. Two of them are the
assertions that could not previously exist and were confirmed RED on the pre-fix `.so`:

- `test_charge_unbalanced_but_element_balanced_is_rejected`: OH → OH⁻ (`electrons=0`), element
  balanced, charge 0 → −1. Pre-fix returned `True`; now `False`.
- `test_sign_inverted_electron_reaction_is_rejected`: the reduction A → A⁻ carried with the
  electron sign flipped (`electrons=+1` instead of `-1`). Pre-fix returned `True`; now `False`.

All five pass after the fix (`5 passed`).

## Hypothesis disposition ("Why this matters")

**Operator's hypothesis:** the dead charge check is why a charge/electron sign inversion recently
reached the mechanism writer undetected, because `_create_reaction` (family.py:1789),
`KineticsLibrary.load` (library.py:567) and `KineticsDepository.load` (depository.py:242) all
gate on `is_balanced()`, and the guard could never fail.

**Disposition: ESTABLISHED, at the mechanism level.** The sign-inverted reduction returned `True`
before the fix (RED test above) and is rejected after it — so the guard demonstrably *could not*
catch a charge/electron sign inversion, and now can. I confirmed the *mechanism* (the guard could
not fire on the class of defect), not the specific historical reaction; I did not re-derive that
exact event. The independent E-balance guard in `electron_balance.py::check_electron_balance`
(untouched, per non-goals) remains a second net at export time — consistent with the campaign's
goal of removing the downstream repair workaround.

## Full-suite fallout — `pytest test/rmgpy/`

`3 failed, 2454 passed, 34 skipped, 14 errors in 219s`. **All 17 failures/errors are confined to
three charged-chemistry *fixture* test files.** Nothing else in the tree regressed (reaction copy
/ pickle / repr, chemkin round-trips, molecule, thermo, solver, rmg, arkane-adjacent — all green).
The named priority suites `reactionTest.py`, `electronPlacementTest.py`, `solver/plasmaTest.py`
are fully green (174 passed, 2 skipped).

### Bucket (c) — test fixtures that assert acceptance of a charge-unbalanced reaction

These fixtures were written to lean on the dead check; several say so verbatim in their comments.
They are *wrong-coverage* and now correctly rejected. They are **not** production chemistry.

**`electronPropagationTest.py` — 10 items** (9 setup-errors + 1 failure). Root cause: the fixtures
stamp `electrons = ±1` onto neutral `CH3 + CH3 <=> C2H6`. This is deliberately non-physical — the
fixture docstring: *"nothing about methyl recombination requires an electron. That is the point."*
Neutral + electron is genuinely charge-unbalanced (charge 0 with a nonzero electron count), so the
live check rejects it at load. The fixtures exercise electron-count *propagation plumbing*, not
chemistry.
- `TestFamilyPropagatesElectronsToItsTrainingDepository::test_family_declaration_is_read`
- `TestFamilyPropagatesElectronsToItsTrainingDepository::test_declared_electrons_reach_the_loaded_training_reactions`
- `TestFamilyPropagatesElectronsToItsTrainingDepository::test_training_reactions_of_an_ordinary_family_are_unchanged`
- `TestFamilyPropagatesElectronsToItsTrainingDepository::test_loaded_training_reactions_have_the_physical_molecularity`
- `TestFamilyPropagatesElectronsToItsTrainingDepository::test_data_borne_count_still_beats_the_family_declaration`
- `TestFamilyPropagatesElectronsToItsTrainingDepository::test_depository_carries_the_count_handed_down_by_its_family`
- `TestReverseTrainingReactionNegatesElectronCount::test_flipped_reverse_training_reaction_negates_electrons` (`C2H6 <=> CH3 + CH3`)
- `TestCreateReactionStoresFamilyForwardOrientation::test_forward_generation_is_family_forward_and_keeps_the_sign`
- `TestCreateReactionStoresFamilyForwardOrientation::test_reverse_generation_is_also_family_forward_and_keeps_the_sign`
- `TestDepositoryDataBorneElectronPrecedence::test_data_borne_count_wins_for_every_charge_transfer_class` (the one FAILED, not ERROR — neutral reaction + `electrons=1` from kinetics)

**`plasmaLocalContextTest.py` — 5 items** (all setup-errors, one bad entry fails the whole class
load). Root cause: entries `Heplus => He` (charge +1 → 0) and `H => proton` (0 → +1) are
electron-implicit. The fixture comment: *"Reaction.is_balanced treats the electron as a conserved
element, so the explicit forms … are both rejected as unbalanced."* Their kinetics
(`BadnellRRArrhenius`, `VoronovEIArrhenius`) set `uses_electron_density` but carry **no signed
electron count**, so nothing can supply the electron to balance them at load. Correctly rejected.
- `test_two_temperature_plasma_entry_loads`
- `test_electron_collision_plasma_entry_loads`
- `test_badnell_rr_arrhenius_entry_loads`
- `test_voronov_ei_arrhenius_entry_loads`
- `test_all_four_entries_loaded_as_library_reactions`

### Bucket (b) — reaction genuinely balanced; the LIBRARY LOADER mis-measures it

**`surfaceChargeTransferBEPLocalContextTest.py` — 2 failures.** `proton + CH2X => CH3X` with
`SurfaceChargeTransferBEP(..., electrons=-1)`. This reaction *is* balanced: H⁺ + CH2X + e⁻ →
CH3X (charge +1, `+= -1` → 0; product 0). The electron count sits in the kinetics, but
`KineticsLibrary.load` (library.py:567) **never folds it into `rxn.electrons` before calling
`is_balanced()`** — unlike `KineticsDepository.load`, which does exactly that at depository.py:234
/ 240 (with the in-code note *"this must happen before `is_balanced()` below inspects the net
charges"*). So the check sees `electrons=0` and correctly rejects the mis-measured reaction.
- `test_surface_charge_transfer_bep_entry_loads`
- `test_get_library_reactions_preserves_electrons`

**This same loader defect is the real production blast radius.** Charged libraries in
`RMG-database-plasma` carry their electron count in the kinetics and load through this path:
- `CO2RR_DFT_Ag111` — 9 `electrons=` entries (e.g. `CO2X + proton <=> CO2HX`).
- `LithiumPrimaryChargedKinetics` — 36 `electrons=` entries (e.g.
  `[Lip] + N=C <=> [Li]N[CH2]` with `ArrheniusChargeTransfer(..., electrons=-1)` — balanced once
  the count is folded in: Li⁺ +1, `+= -1` → 0).
These are exercised by `@pytest.mark.database` tests (excluded from the unit run above) and by any
full RMG database load. With the check live and the loader unfixed, **every one of them breaks**.

### Bucket (a) — genuinely charge-unbalanced production chemistry that should always have been rejected

**None found** in the unit suite. The electron_propagation fixtures are non-physical by design
(bucket c); the plasma-local-context entries are wrong-coverage fixtures (bucket c); the
production charged libraries are all balanced-but-mis-measured (bucket b). A database-marked and
regression sweep would be needed to close bucket (a) across the whole shipped database.

## Recommended staging (Remaining Work)

The comparison is landed and correct and **must not be narrowed**. The remaining work is a
separate, reviewable change because it touches production database loading and fixture design:

1. **Measurement fix (bucket b), in scope, well-defined:** mirror `KineticsDepository.load`'s
   electron propagation into `KineticsLibrary.load` — read the kinetics' `electrons` count (and,
   where applicable, the family declaration) onto `rxn.electrons` *before* `library.py:567`'s
   `is_balanced()`. This rescues `surface-charge-transfer-bep` and every production charged
   library (CO2RR, Lithium*). It is not narrowing — it makes the loader measure charge the same
   correct way the depository already does. Needs its own review: it changes how all charged
   libraries load, and may itself surface genuine bucket-(a) reactions where a library's declared
   count does not actually balance.
2. **Fixture reworks (bucket c):** the `electron_propagation_data` and `plasma-local-context`
   fixtures must be redesigned to test their plumbing with reactions that actually balance (charged
   species, or a representation carrying a real signed electron count) rather than neutral
   stand-ins / electron-implicit forms that leaned on the dead check. For the plasma entries this
   is a design question — `BadnellRRArrhenius`/`VoronovEIArrhenius` carry no signed electron count
   today, so representing electron-implicit recombination/ionization under a live charge check
   needs a decision, not a mechanical edit.
3. **Database + regression sweep** to close bucket (a) across the shipped `RMG-database-plasma`.
