# I-208 — let a conserved-electron reaction resolve its electron placement

**Worktree:** `/home/alon/Code/RMG-Py-i208-gate`, branch `i208-gate`, cut from `plasma` (0 behind).
**Database (read-only):** `/home/alon/Code/RMG-database-plasma` @ `plasma`.
**Nothing pushed, nothing merged.** Adversarial review expected before this lands.

## The one-line change

`PlasmaReactor._resolve_electron_placements` selected reactions for the electron-placement resolver
with `if getattr(rxn, 'electrons', 0):` alone (`rmgpy/solver/plasma.pyx:379`). That gate closes on a
net of zero, so a **conserved**-electron reaction (`AB + e- -> A + B + e-`, one electron in, one out)
never reached the resolver and kept its heavy-only molecularity while its coefficient is second
order — a first-order evaluation of a second-order rate, wrong by one factor of the electron density.
The gate now also admits a net-zero reaction whose owner declares a placement and which carries no
explicit electron yet. Everything else still passes through untouched.

## Questions 1–4 (measured)

### 1. Reproduce the defect — the resolver is not consulted for a net-zero carrier

`test/rmgpy/electronPlacementTest.py::TestNetZeroElectronPlacementGate` builds a canonical
conserved-electron reaction `AB -> A + B` (`electrons=0`, no explicit electron, second-order
`TwoTemperaturePlasma` coefficient `cm^3/(mol*s)`) whose owner is monkeypatched into
`FAMILY_ELECTRON_PLACEMENT` as `(1, 1)` — the same technique `i113IonisationPlacementTest.py` uses to
exercise the code path without shipping the data. Before the fix, run against the OLD `.so`:

```
3 failed, 3 passed
  FAILED …test_net_zero_carrier_reaches_the_resolver        # returned by identity: gate skipped it
  FAILED …test_resolved_order_matches_the_coefficient_units
  FAILED …test_resolver_is_actually_consulted   assert 0 == 1  (spy never called)
```

The spy on `electron_placement.resolve_electron_placement` records **zero** calls: the reaction is
passed through with one reactant against an order-2 (`cm^3/(mol*s)`) coefficient. That is the defect —
first order applied to a second-order rate.

### 2. Every gate, not just `:379`

`electrons` non-zero-ness decides placement in exactly one place — the `:379` selection gate. The
sibling references are not gates of the same kind:

- `plasma.pyx:676` (`_validate_reactions`) *rejects* a metadata-only electron count; it never selects
  for placement. After resolution a net-zero view carries `electrons=0` and an explicit electron, so
  it passes this guard. Unchanged.
- `plasma.pyx:692`/`:739` reject a plasma-flagged rate law with no explicit incident electron. These
  are the **reactor backstop**: for a net-zero family with `TwoTemperaturePlasma`
  (`uses_electron_temperature=True`, `uses_electron_density=False`), `:739` would have loudly refused
  the reaction in a full reactor run — so the "silent first-order rate" is silent only at the
  database/generation layer, not through `PlasmaReactor`. The proper fix is to *place* the electron
  (this ticket), not merely to reject; the backstop stays as defence in depth.
- The resolver itself (`electron_placement.py`) needs **no change**. Tracing `(1,1)`, `electrons=0`:
  step 5 `expected_net = product_count - reactant_count = 0` matches; step 10 E-balance holds with
  one electron each side; step 11 order `= len(view.reactants) = 2` matches the `cm^3/(mol*s)`
  coefficient. The export path (`expand_electrons`, Chemkin, Cantera) already reads
  `FAMILY_ELECTRON_PLACEMENT` per `Reaction.family` (I-126) and needs no gate change — it places by
  declaration, not by net.

### 3. What the gate protects, and the correct predicate

The old gate skipped the resolver for ordinary neutral reactions — the overwhelming majority, which
have no electron anywhere and no declaration. Running the resolver on them would raise
`ElectronPlacementError` by name for every one. The predicate must keep skipping those. Chosen
predicate (`_needs_electron_placement`):

```
electrons != 0                                    # metadata electron: historical case, unchanged
  OR (owner in FAMILY_ELECTRON_PLACEMENT           # a placement is declared for this owner …
      AND no explicit electron already present)    # … and it is not already in reactor form
```

- The metadata clause is kept verbatim, so an undeclared metadata-electron reaction still fails by
  name inside the resolver, and a double representation is still caught at resolver step 2 — no
  current behaviour changes there.
- The declaration clause is the widening. It uses the same `get_placement_owner` the resolver keys
  on, so gate and resolver cannot disagree about ownership. The `no explicit electron` conjunct keeps
  an already-reactor-form reaction (explicit electron, `electrons=0`) passing through by identity
  rather than hitting the resolver's "nothing to resolve" refusal.

Not chosen: the bare "an electron appears among its species" candidate, because routing an
explicit-electron reaction to the resolver makes step 4 raise "already carries explicit electron;
nothing to resolve" — a false failure for a legitimately reactor-form reaction.

### 4. Blast radius — who newly reaches the resolver

Newly-reaching = `electrons == 0` **and** declared owner **and** no explicit electron. Every owner
currently in `FAMILY_ELECTRON_PLACEMENT` declares a net != 0 pair (`(1,0)`, `(1,2)`, `(0,1)`), so a
well-formed reaction of those owners already carries `electrons != 0` and was already routed. The
only net-zero declared reaction reachable **today** is a malformed one (declared owner, zero net),
which the resolver's step 5 refuses by name — the correct loud failure. So the newly-reaching set
against the shipped registry is empty for well-formed input; confirmed empirically by all standing
suites green and the argon deck unchanged. When the separate data ticket adds a `(1,1)` entry, its
reactions newly reach and resolve — the intended outcome.

## Docstring correction

`rmgpy/electron_placement.py`'s module docstring attributed
`Plasma_Electron_Impact_Dissociation`'s absence to a **spectator** the family model cannot express.
That is stale for the net-zero shape: the recipe (`BREAK_BOND *1-*2`, `GAIN_RADICAL *1`,
`GAIN_RADICAL *2`) touches both heavy centres, so no heavy spectator arises, and the electron is
absent from the template by convention. The real blocker was the `:379` gate — now lifted. The
docstring now says so, and records that the family is still out of the table only because no `(1,1)`
declaration has been added yet (a separate data ticket). `Plasma_Collisional_Ionization`'s spectator
note is left intact (out of scope; not re-measured here).

## Files changed

- `rmgpy/solver/plasma.pyx` — `_resolve_electron_placements` calls the new `_needs_electron_placement`
  predicate; docstrings updated.
- `rmgpy/electron_placement.py` — module docstring corrected (spectator → gate), two passages.
- `test/rmgpy/electronPlacementTest.py` — `TestNetZeroElectronPlacementGate` (6 tests) + one import.

## Non-goals honoured

No placement declaration added to `FAMILY_ELECTRON_PLACEMENT`; the two-number map's meaning,
validation rule, and existing entries untouched; RMG-database not modified; nothing else repaired;
nothing pushed or merged.

## Reported, not fixed

`rmgpy/molecule/`, `rmgpy/kinetics/`, `rmgpy/data/` are gated and were not entered. The reactor
backstop asymmetry noted in Q2 (a net-zero *undeclared* plasma family is refused at `:739` only when
its kinetics set `uses_electron_temperature`/`uses_electron_density`; a plain-Arrhenius net-zero
electron family would still be evaluated first-order silently) is a real residual gap, but it is
outside this gate's remit — it needs either a declaration for such a family or a broader
`_validate_reactions` guard, both separate tickets.
