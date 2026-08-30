# I-169 — what a modeller types to ask for a charge-neutral plasma

**Status: PROPOSAL, not a ruling.** The spelling of this keyword is a public input-file API. Once a
deck can type it, this project supports it indefinitely and cannot quietly rename it. It was
surfaced as `NEEDS-INPUT` before implementation; option A is implemented on the branch so the rest
of the ticket could be verified, and is labelled a proposal throughout.

Written against `plasma@4f18dc389`, worktree `/home/alon/Code/RMG-Py-i169-charge-neutral`.

## The four questions every candidate must answer

1. What does the modeller type?
2. What happens when a cation is already declared?
3. What happens when several cations could balance the charge?
4. What happens to a deck that says nothing?

## A — `chargeBalanceSpecies='Ar+'` · **recommended**

```python
plasmaReactor(
    temperature=(298.15, 'K'),
    pressure=(5, 'torr'),
    electronTemperature=(34813.5, 'K'),
    electronDensity=(1e16, 'm^-3'),
    initialMoleFractions={'Ar': 1.0},
    chargeBalanceSpecies='Ar+',          # <- the whole change
    terminationTime=(1e-3, 's'),
)
```

1. **Types**: a keyword on `plasmaReactor(...)` naming the label of a declared species that absorbs
   whatever charge the rest of the composition leaves over. RMG computes its mole fraction.
2. **Cation already declared**: declaring it is *required* — the label must resolve to a species the
   deck declared. What is refused, with a named `InputError`, is naming it **and** typing its mole
   fraction in `initialMoleFractions`: those are two sources of truth for the same number and may
   disagree. This is the exact shape `electronDensity` already refuses for the electron
   (`input.py:631-637`), so the rule is one modellers have already met.
3. **Several cations**: no ambiguity can arise, because the deck names exactly one. Any other
   charged species may carry explicit fractions; they enter the balance, and the named species takes
   the remainder. This is the decisive advantage over B.
4. **Deck says nothing**: composition identical to today's, to the last bit. It additionally gets
   the new net-charge warning if it is non-neutral.

**Why this spelling.** The input DSL already has `balanceSpecies=<label>` on `simpleReactor` and
`constantVIdealGasReactor`, meaning precisely "the species whose mole fraction RMG computes instead
of the modeller typing it". `chargeBalanceSpecies` is that idea, qualified by *what* is being
balanced. Same camelCase, same label-valued shape, no new mental model.

**Cost.** One more keyword, and the modeller must decide which ion carries the compensating charge —
which is a modelling decision they should be making explicitly anyway, not one RMG should guess.

## B — `chargeNeutral=True`

1. **Types**: a boolean.
2. **Cation already declared**: fine if there is exactly one.
3. **Several cations**: **this is the killer.** RMG must *choose*, and every available rule —
   alphabetical, first-declared, lowest ionisation energy, largest existing fraction — is arbitrary.
   Worse, it is silently arbitrary: a deck that works correctly today changes its physics the day it
   gains a second cation, with no message. And with *zero* cations declared the flag has no answer
   at all: it must either raise (so the flag cannot stand alone, and the modeller is back to naming
   a species) or invent a species, which is never acceptable.
4. **Deck says nothing**: unchanged.

**Cost.** Cheapest to type; a trap that fails silently as decks grow. A flag that is only correct in
the one-cation case should not be the public API for a general one.

## C — sentinel in the composition, `initialMoleFractions={'Ar': 1.0, 'Ar+': 'balance'}`

1. **Types**: a magic string where a float belongs.
2/3. As A — the sentinel names the species, so no ambiguity.
4. Unchanged.

**Cost.** Type-punning a float-valued dict. Every consumer of `initialMoleFractions` must learn the
sentinel: the finiteness and sign guards (`input.py:585-596`), the ranged-condition detector (which
already gives list values a meaning), and `save_input_file`'s round-trip. A typo — `'balanced'` —
becomes an unnamed `float()` crash instead of a named `InputError`. Reads nicely; costs more than it
reads.

## D — `neutralizeWith='Ar+'`

Semantics identical to A. **Cost:** "neutralize" is overloaded in chemistry (acid–base), and there
is no precedent for the shape in this DSL. `chargeBalanceSpecies` says *charge* out loud.

## E — `chargeCarrier='Ar+'`

**Cost, and it disqualifies the name:** in plasma physics "charge carrier" means every mobile
charged species, electrons included — not the one species that closes the balance. Actively
misleading to the exact audience that will read it.

## The sub-decision, also surfaced

Does the keyword work **only** alongside `electronDensity`, or also when the electron fraction is
typed explicitly in `initialMoleFractions`?

**Recommendation: both.** In the explicit-electron branch the arithmetic is four lines (the
balancing fraction is just `-Σ x_i z_i / z_balance` before normalisation), and a public keyword that
silently does nothing in one of the two branches is a defect being shipped on purpose. Both branches
are implemented.

## Not part of this proposal

The **warning threshold** is an internal numeric constant, not a public API, so it is settled in the
code and documented there rather than escalated: `PLASMA_NET_CHARGE_ATOL = 1e-12` (absolute floor,
above the `N·eps` double-precision accumulation bound for models up to ~4500 species) and
`PLASMA_NET_CHARGE_RTOL = 1e-6` (relative to `Σ |x_i z_i|`, the total charge magnitude the net is
measured against). The check warns when `|net| > max(ATOL, RTOL · Σ|x_i z_i|)`. Both terms are
needed: a purely absolute test cannot see a weakly-ionised deck that is 100% imbalanced, and a
purely relative test cries wolf when the charged fraction is near the roundoff floor.
