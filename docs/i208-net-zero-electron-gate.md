# I-208 — let a conserved-electron reaction resolve its electron placement

**Worktree:** `/home/alon/Code/RMG-Py-i208-gate`, branch `i208-gate`, rebased onto `plasma`
@ `8751394fe` (after I-206 landed the two family declarations).
**Database (read-only):** `/home/alon/Code/RMG-database-plasma` @ `plasma`.
**Nothing pushed, nothing merged.** Manager confirmed the gate merges (D-048, round 16); this revision
answers the round's four follow-ups.

## The one-line change

`PlasmaReactor._resolve_electron_placements` selected reactions for the electron-placement resolver
with `if getattr(rxn, 'electrons', 0):` alone (`rmgpy/solver/plasma.pyx:379`). That gate closes on a
net of zero, so a **conserved**-electron reaction (`AB + e- -> A + B + e-`, one electron in, one out)
never reached the resolver and kept its heavy-only molecularity while its coefficient is second
order — a first-order evaluation of a second-order rate, wrong by one factor of the electron density.
The gate now also admits a reaction whose owner declares a placement, whose `electrons` scalar is 0,
and which carries no explicit electron yet. The net-zero conserved-electron reaction is the
**motivation** for that clause, not its boundary: because the clause keys on `electrons == 0` for any
declared owner, it also admits a declared **nonzero** owner (e.g. `Plasma_Electron_Attachment` at
`(1, 0)`) arriving with a defaulted or corrupt zero `electrons` scalar and no explicit electron. Such
a reaction now hard-fails **by name** at the resolver's net-mismatch check instead of passing through
silently as the old one-line gate let it — that loud failure is intentional. Everything outside this
clause (an ordinary neutral with no electron and no declaration; a reaction already carrying its
electron explicitly) still passes through untouched.

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
  coefficient. **The export path is a separate concern with its own defect for the net-zero shape**,
  corrected below in the addendum (round 16): `expand_electrons` early-returns at `electrons == 0`
  *before* it reads the declaration, so it does NOT place by declaration for a net-zero reaction.
  I-208 does not change the export path; the round established that the export path also needs work,
  as a HARD PREREQUISITE for the data ticket — see the addendum.

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
That is stale for the net-zero shape. The real blocker was the `:379` gate — now lifted. The docstring
now says so, and records that the family is still out of the table only because no `(1,1)` declaration
has been added yet (a separate data ticket). `Plasma_Collisional_Ionization`'s spectator note is left
intact (out of scope; not re-measured here).

**Recipe provenance (round-16 REPORT-4 — this reason has inverted twice, so the citation is exact).**
The "touches both heavy centres, no heavy spectator" claim is read directly from the family's own
recipe. It is **not** in RMG-Py, and **not** in the database's active `input/` tree — the family is
staged out — which is why a reviewer standing in RMG-Py could not check it. It lives in the
**RMG-database `plasma` branch**, at
`docs/plasma_electron_impact_dissociation_staged/groups.py`. Verbatim:

```python
template(reactants=["AB"], products=["A", "B"], ownReverse=False)   # line 18
...
electrons = 0                                                        # line 24
recipe(actions=[                                                     # lines 26–30
    ['BREAK_BOND', '*1', 1, '*2'],
    ['GAIN_RADICAL', '*1', 1],
    ['GAIN_RADICAL', '*2', 1],
])
```

Both labelled centres `*1` and `*2` are members of the single reactant group `AB`, and every action
operates on one or both of them (`BREAK_BOND *1-*2`, `GAIN_RADICAL *1`, `GAIN_RADICAL *2`). There is
no third labelled atom, so no heavy participant is a spectator; and `electrons = 0` with no electron
in the template confirms the net-zero, template-implicit electron. The claim is therefore verified at
its source, and the docstring now carries this path so the next reader can check it without leaving
the code.

**One artifact to NOT copy from (round-17 ITEM 4).** A second staged copy exists at
`docs/i154-carry-chemistry/held-back/Plasma_Electron_Impact_Dissociation/groups.py`, and it is **not
equivalent** — its template is `template(reactants=["AB", "e-"], products=["A", "B", "e-"])` with the
electron declared as an explicit template participant (entry label `"e-"`, `1 *3 e u1 p0 c-1`). A
reaction generated from *that* form carries an explicit electron and therefore **bypasses this gate
entirely** (the gate passes explicit-electron reactions through as already-reactor-form). The
`plasma_electron_impact_dissociation_staged` copy cited above — electron template-implicit,
`electrons = 0` — is the one the net-zero gate is designed for. Anyone adding the `(1,1)` entry must
copy from the template-implicit copy, not the explicit-electron one.

## Files changed

- `rmgpy/solver/plasma.pyx` — `_resolve_electron_placements` calls the new `_needs_electron_placement`
  predicate, which uses the resolver's defensive `_is_electron` helper (round-16 MUST-FIX 2);
  docstrings updated.
- `rmgpy/electron_placement.py` — module docstring corrected (spectator → gate) with the recipe
  citation; two passages, rebased to co-exist with I-206's two-family note.
- `test/rmgpy/electronPlacementTest.py` — `TestNetZeroElectronPlacementGate` (6 tests) + one import;
  reassembled after I-206's tests to keep both suites.

## Non-goals honoured

No placement declaration added to `FAMILY_ELECTRON_PLACEMENT`; the two-number map's meaning,
validation rule, and existing entries untouched; RMG-database not modified; nothing else repaired;
nothing pushed or merged.

## Addendum — the third call site (adversarial review of I-206)

A review flagged a third consumer of the placement declarations,
`rmgpy/electron_balance.py:284 get_electron_placement_counts()`. Verified against the code:

**(a) Site 3 is genuinely ungated, and a `(1,1)` resolves there for a net-zero reaction.**
`get_electron_placement_counts` calls `get_placement_declaration(reaction)` at line 284 with **no**
`electrons != 0` guard. For a net-zero reaction whose owner declares `(1,1)`: `declaration=(1,1)`,
`electrons=0`, and line 287 `product_count - reactant_count == electrons` → `1 - 1 == 0` → returns
`(1, 1)`. Confirmed. So the blanket premise "a net-zero declaration could never be consulted anywhere"
is **false** — site 3 consults it today. (The I-208 code docstring's narrower statement is scoped to
`resolve_electron_placement`/`_resolve_electron_placements` — "the resolver" — and remains accurate:
that consumer, site 1, was the one gated out.)

The three sites, and what each answers for a net-zero reaction:

| # | Site | Gate on net-zero | Net-zero `(1,1)` declared |
|---|------|------------------|---------------------------|
| 1 | `plasma.pyx` `_resolve_electron_placements` (reactor) | **was** `electrons != 0`; I-208 lifts it | now resolves → places `(1,1)` |
| 2 | `electron_balance.py` `expand_electrons` (export) | `if electrons == 0: return` early (line 351) | **early-returns, places nothing** |
| 3 | `electron_balance.py` `get_electron_placement_counts` (identity) | **none** | consults declaration → `(1,1)` |

**(b) Consequence of the current state (no net-zero owner declared) — a distinct, latent defect.**
`are_identical_species_references` (`rmgpy/rmg/model.py:2356`) compares heavy-species references (in
either direction) **and** the two `get_electron_placement_counts` pairs — and nothing else (not
family, not kinetics). With no net-zero declaration, a catalytic-electron reaction falls through site
3 to `(0, 0)`. A thermal `AB => A + B` and a plasma net-zero `AB => A + B` then present identical
heavy references **and** identical `(0,0)` counts, so the predicate declares them the same reaction
and the model silently keeps whichever was offered first.

- **Reachable?** The predicate has no family/kinetics discriminator, so the collapse is reachable *in
  principle*. It is **not** reachable with shipped data: it needs a plasma net-zero dissociation
  family/library to coexist with a thermal reaction over the same heavy species, and no such family
  is in the database (i204 staged `Plasma_Electron_Impact_Dissociation` out of `input/`).
- **Distinct from I-208?** Yes — this is an **identity collapse at the model builder**, separate from
  I-208's rate-**order** symptom at the reactor. Same root cause, though: the missing `(1,1)`
  declaration. Adding it makes site 3 return `(1,1) ≠ (0,0)`, which distinguishes the two reactions —
  so the separate data ticket closes this at the same stroke it closes the order defect. Named here;
  **not fixed** (touches `rmgpy/rmg/model.py`, and fixing it without the declaration is untestable).

**(c) Does `_needs_electron_placement` make sites 1/2/3 agree? No — and I am not widening scope to
force it.** The predicate is a site-1-only change.

- **Today (no net-zero declaration): all three agree.** A well-formed net-zero reaction is "no
  placement" at every site — site 1 passes it through (owner undeclared, so the new declaration
  conjunct is false), site 2 early-returns, site 3 returns `(0,0)`. My change does not disturb this:
  it routes a net-zero reaction only when its owner is declared, of which there are none, plus a
  *malformed* declared-owner-with-zero-net reaction, which the resolver then refuses by name.
- **When the data ticket adds `(1,1)`: sites 1 and 3 agree (both place `(1,1)`); site 2 does not.**
  `expand_electrons` retains its `electrons == 0` early return, so it exports `AB => A + B` with no
  electron and a bimolecular coefficient.

**CORRECTION (round 16) — my earlier claim that this fails loudly was WRONG, and I verified the
correction against the code myself.** I had written that the export path pairs `expand_electrons` with
`check_electron_reactant_order`, which would raise `MechanismWriterError` on the order mismatch. That
is false for exactly the kinetics a net-zero dissociation family uses, and two export paths never call
the guard at all:

1. **The order guard goes silent for generic Arrhenius.** `check_electron_reactant_order`
   (`rmgpy/electron_balance.py:503`) has two raise paths. The first needs `kinetics.uses_electron_density`
   (only `BadnellRRArrhenius`/`VoronovEIArrhenius` set it). The second computes
   `required_order = get_plasma_rate_order(kinetics)` and then `if required_order is None or
   required_order == len(reactants): return`. `get_plasma_rate_order` (`:479`) returns a number only
   for `TwoTemperaturePlasma`, `ElectronCollisionPlasma`, `BadnellRRArrhenius`, `VoronovEIArrhenius`;
   for a generic `Arrhenius`/`ArrheniusEP` it returns **`None`**, so the guard **returns silently**.
   My test used `TwoTemperaturePlasma` (order-recognised), which is precisely why the guard fired in
   my head — I generalised "loud" from the one kinetics class that happens to be covered. The two
   data-less families I-206 shipped carry **`ArrheniusEP`** rules with `m^3/(mol*s)` units, and a
   future `(1,1)` dissociation family is the same shape — the exact blind spot.
2. **`Reaction.to_cantera()` never calls the guard.** `rmgpy/reaction.py:319` folds electrons via
   `expand_electrons` and builds the in-memory Cantera reaction, with no `check_electron_reactant_order`
   anywhere in the method (verified: the only calls are in `chemkin.pyx` and `yaml_cantera2.py`).
3. **`CanteraWriter1` walks past its front gate.** `yaml_cantera1.py:138` refuses a mechanism only when
   `next((rxn for rxn in rxns if getattr(rxn, 'electrons', 0)), None)` is truthy. A conserved-electron
   reaction has `electrons == 0`, so it is not selected and the writer emits `AB => A + B` via
   `_build_equation_string` from the raw participant lists — electron omitted, silently.

So a future net-zero family would export `AB => A + B` with a bimolecular coefficient and no electron,
**silently** — the exact defect class I-208 exists to prevent, one layer out.

**HARD PREREQUISITE for the `(1,1)` data ticket (recorded so it is not missed).** Adding the
declaration alone is a silent wrong export. Whoever adds it must ALSO, in the same change:
- **(P1)** lift `expand_electrons`'s `electrons == 0` early return so it consults the declaration for a
  net-zero owner (placing the electron on both sides in the declared numbers); and
- **(P2)** extend `check_electron_reactant_order` to recognise a generic `Arrhenius`/`ArrheniusEP`
  whose `A.units` imply an order higher than the exported reactant count — i.e. read the order from the
  units, not only from the four plasma classes — so the guard stops being blind to the shape the data
  ticket will actually use. `to_cantera()` and `CanteraWriter1` should route through the same guard or
  be documented as unsupported for plasma mechanisms.

These are inert today (no net-zero declaration exists) and cannot be tested end to end without the
declaration, which is why they belong with the data ticket and not here — but they are its
preconditions, not optional follow-ups.

**Ruling — the narrow site-1 fix stands, and the deferral now rests on the correct reason.** Site 3 is
not broken; it is the reference-correct behaviour (ungated, already honours net-zero declarations), and
changing it would be the regression. Site 2's early return, its silent export blind spots, and the
site-(b) identity collapse are all gated on the same *missing declaration*, are inert with shipped
data, and belong with the data ticket that adds `(1,1)` — carried above as hard prerequisites P1/P2 so
that ticket cannot land a silent wrong export. I-208 delivers one necessary precondition (site 1) and
does not claim more; my earlier "fails loud" safety claim was the round's real finding and it was
wrong.

## Round-16 MUST-FIX 2 — the gate uses the resolver's defensive electron test

`_needs_electron_placement` called `spc.is_electron()` behind a `callable(getattr(...))` guard, which
catches a *missing* method but not an `IndexError` from `spc.molecule[0]` on a malformed participant —
so a declared-owner net-zero reaction with a malformed species failed **at the gate** with a raw
`IndexError` instead of reaching the resolver and failing with the named `ElectronPlacementError`. Now
the gate calls the resolver's own `electron_placement._is_electron` (`rmgpy/electron_placement.py:250`),
which catches `(AttributeError, IndexError)`. Gate and resolver now agree on what an electron is *and*
fail the same way. Folded into the gate commit.

## Round-16 REPORT-3 — two more shapes the round raised (decided, not fixed)

**(a) Partial explicit form `AB -> A + B + e-` with `electrons == 0` and a declared owner.** The gate
treats *any* explicit electron as already-reactor-form and passes it through, so this one-sided form is
neither resolved nor refused. **Decision: no gate change.** The gate cannot tell a well-formed
reactor-form reaction (balanced explicit electrons, the reachable case — a resolver view fed back)
from a one-sided partial (unreachable from the resolver or normal generation) without re-implementing
the resolver's balance logic inside the gate — which is precisely the multi-site-disagreement smell
this round is about. For the reachable case (balanced) pass-through is correct; for the partial case,
if its kinetics is plasma-flagged the reactor still refuses it at `_validate_reactions:739`
(`n_electron_reactants == 0`), and if it is generic Arrhenius it is the same silent-order gap already
recorded under the export prerequisites — not a new one. Optimising the gate for the unreachable case
would regress the reachable one.

**(b) Species-less declared net-zero reaction (`reactants=[]`, `products=[]`).** Under the predicate
it is routed (declared owner, no explicit electron), and the resolver builds placement into `e- => e-`,
which reactor validation then accepts — malformed input made runnable. **Decision: no gate change;
reported as a resolver-hardening gap.** The object that *builds* `e- => e-` from empty lists is the
resolver, not the gate; the right refusal ("a placement that leaves no heavy participant is
malformed") belongs there (or at reaction construction), and adding it is resolver-scope, outside this
ticket's remit and its non-goals. A degenerate empty-participant reaction is not produced by RMG
generation; it is hand-constructed. Putting an emptiness check in the routing gate would, again,
duplicate validation across sites. Flagged for the resolver-hardening owner.

## Round-17 tightening (D-048)

- **ITEM 1 — the gate is wider than net-zero, kept and documented.** `_needs_electron_placement`
  admits `electrons == 0` for **any** owner in `FAMILY_ELECTRON_PLACEMENT`, not only owners whose
  declared net is zero. So a declared **nonzero** owner (e.g. `Plasma_Electron_Attachment` at `(1,0)`)
  that arrives with `electrons = 0` and no explicit electron — a corrupt reaction — now routes to the
  resolver and hard-fails by name at the net-mismatch check (declared net `-1` ≠ metadata `0`), where
  the old one-line gate passed it through **silently**. That is a strict improvement (loud failure of
  corrupt input) and is now stated as intentional in the `_needs_electron_placement` docstring. The
  width is deliberately **not** narrowed with a `product_count == reactant_count` condition — that
  would restore the silent pass-through for exactly the corrupt input it should catch. New tripwire
  test: `test_declared_nonzero_owner_with_zero_metadata_fails_by_name`.
- **ITEM 2 — the defensive-gate test now pins the follow-through.**
  `test_gate_electron_test_is_defensive_like_the_resolver` no longer asserts only that the gate returns
  `True`; it now calls `_resolve_electron_placements` and asserts the malformed participant reaches the
  resolver and raises `ElectronPlacementError` (the E-balance step wraps the `IndexError`,
  `electron_placement.py:646–650`) — the named failure the gate defers to, actually observed.
- **ITEM 3 — "silent" narrowed in code/test comments.** The `_resolve_electron_placements` docstring
  and the test comment block now say explicitly that on the full reactor path with plasma-flagged
  kinetics the wrong order is caught at `_validate_reactions` (which RAISES), and that the genuinely
  silent form is the export path with generic Arrhenius-like kinetics — matching this report rather
  than implying the reactor mis-evaluates silently.
- **ITEM 4 — the wrong-to-copy staged artifact is called out.** See the Docstring-correction section:
  the second staged copy
  (`docs/i154-carry-chemistry/held-back/Plasma_Electron_Impact_Dissociation/groups.py`) declares the
  electron as an **explicit** template participant, so a reaction from it would bypass this gate
  entirely; the `(1,1)` data ticket must copy from the template-implicit
  `plasma_electron_impact_dissociation_staged` copy.

## Pre-existing defect surfaced (round-17, NOT fixed here — reported)

**Seed-mechanism reload drops a reaction's metadata electron count.** `make_seed_mech`
(`rmgpy/rmg/main.py:2051,2072`) writes each seed entry's `label = reaction.to_labeled_str()`, and
`to_labeled_str` (`rmgpy/reaction.py:199`) builds the equation from participant labels only — no
electron folding. On the seed write path, `save_entry` (`rmgpy/data/kinetics/common.py:75-90`)
serializes the reaction item as an explicit field whitelist — `degeneracy`, `duplicate`,
`reversible`, `allow_pdep_route`, `elementary_high_p`, `allow_max_rate_violation` — and `electrons` is
not among them; the item is never `repr`'d on this path, so its metadata count is not written. (This
is the specific path that loses it, not a blanket absence: `LibraryReaction.__repr__`
(`rmgpy/data/kinetics/library.py:145`) *does* emit `electrons=` when nonzero, so a code path that
serialized `repr(entry.item)` would preserve it — the seed path is not that path.) So a seed reload,
reconstructing the equation from the electron-free `label`, rebuilds the
reaction as `Reaction(electrons=0)` and recovers the electron order **only** when `entry.data` (the
kinetics) carries its own electron information: `VoronovEIArrhenius`/`BadnellRRArrhenius` do (they set
`uses_electron_density` etc.), so they survive the round trip; a **generic family estimate** that
carried its electron as a metadata `electrons` scalar does not — it reloads as neutral, first order.
This is independent of the I-208 gate and pre-exists it; surfaced while tracing the placement
consumers. Owner: whoever holds seed-mechanism serialization. Not touched here (main.py / data/ are
outside this ticket).

## Reported, not fixed

`rmgpy/molecule/`, `rmgpy/kinetics/`, `rmgpy/data/` are gated and were not entered. The reactor
backstop asymmetry noted in Q2 (a net-zero *undeclared* plasma family is refused at `:739` only when
its kinetics set `uses_electron_temperature`/`uses_electron_density`; a plain-Arrhenius net-zero
electron family would still be evaluated first-order silently) is a real residual gap, but it is
outside this gate's remit — it needs either a declaration for such a family or a broader
`_validate_reactions` guard, both separate tickets.
