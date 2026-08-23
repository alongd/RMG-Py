# I-085 — Walls A and B cleared (durable findings)

**Date:** 2026-08-23 · **Branch:** `i085-compose-walls` · **Base:** `6d3c03823`
**Worktree:** `/home/alon/Code/RMG-Py-i085-compose` · **Database:** `/home/alon/Code/RMG-database-plasma/input`

Scope: the first two of the four M7-preflight walls between a neutral lithium-bearing plasma input
and a written mechanism. Walls C and D deliberately untouched. Reproduction decks and verbatim
tracebacks under `input.py`, `input_secondary.py`, `probes/`, `repro-evidence/`.

## Wall A — the electron's spin multiplicity (fixed)

**Decision: the `u1` declaration is the non-canonical side; the fix makes the electron thermo
lookup multiplicity-agnostic and refuses to hand the electron to group additivity / RDKit.**

Evidence:
- `Atom.is_electron()` is `self.element.number == -1` (`rmgpy/molecule/molecule.py:363`) — deliberately
  multiplicity-agnostic. An electron is an electron whether declared `u0` or `u1`.
- RMG's canonical electron representation is `u0`: the SMILES-lookup dict in
  `rmgpy/molecule/translator.py` maps `'e'` → `1 e u0 p0 c-1`, and the plasma database thermo entry
  (`thermo/libraries/electrocatThermo.py`, label `"electron"`) is `1 e u0 p0 c-1`. Both agree on `u0`.
- A free electron physically carries spin ½, but RMG models it as a *charge-carrier pseudo-species*
  with a single canonical thermo entry (H298 = 0), not a spin-bearing radical — which is exactly why
  `is_electron()` ignores multiplicity. The thermo library match, however, keyed on
  `molecule.is_isomorphic(entry.item)` (`rmgpy/data/thermo.py:2057`), which distinguishes `u0` from
  `u1`, so a `u1` electron missed the `u0` entry, fell through to group additivity, and crashed inside
  RDKit (`RuntimeError: ... Element 'e' not found`).

Fix (all in `rmgpy/data/thermo.py`, pure Python — no rebuild):
1. `get_thermo_data_from_library`: when the query species is an electron, match it against any
   electron library entry regardless of unpaired-electron count. This makes `u0` and `u1` both
   resolve to the canonical entry — the lookup succeeds, per the prompt's "make the lookup find the
   right entry, do not make the RDKit fallback succeed."
2. `get_thermo_data`: intercept the electron before the group-additivity fallback. If a library
   entry was found, return it; otherwise raise a clear `DatabaseError` telling the user to load an
   electron thermo library — the electron must never be handed to RDKit.
3. `rmgpy/rmg/input.py:614`: the `electronDensity`-missing-electron error message suggested declaring
   the electron as `u1`; corrected to the canonical `u0`.

## Wall B — the unguarded `Trange` (fixed)

**Decision: give `PlasmaReactor` the `Trange` attribute (conform to the interface), not guard the
caller.**

Evidence:
- `Trange` is a declared public attribute on *every* other concrete reactor — `simple.pyx:110`,
  `liquid.pyx:61`, `surface.pyx:59` each `cdef public list Trange`, assigned only when `T` is a list
  (ranged) and left `None` (its C default) for scalar-`T` conditions. `PlasmaReactor` was the only
  concrete reactor that omitted the declaration.
- `RMG.execute()` (`main.py:879-880`) reads `x.Trange[0] if x.Trange else x.T` on every reaction
  system to compute `Tmin/Tmax`. The `None` branch already falls back to the scalar `x.T`.
- `PlasmaReactor` definitionally forbids ranged T/P/Te (`plasma.pyx:105-108` raises
  `PlasmaStateError` on a list), so its `Trange` is permanently and truthfully `None`. The scalar
  gas temperature `self.T` is the meaningful `Tmin/Tmax` contribution — not a placeholder lie, but
  the identical value a scalar `SimpleReactor` carries.

Fix: declare `cdef public list Trange` on `PlasmaReactor` (`rmgpy/solver/plasma.pyx`), left
unassigned → `None`. Requires a Cython rebuild (`make build`).

## Deck-config gap surfaced (reported, not one of the four walls)

The committed preflight **primary** deck (`input.py`) loaded thermo libraries
`LithiumPrimaryThermo, LithiumAdditionalThermo, primaryThermoLibrary` — none of which contain an
electron entry — and did **not** load `electrocatThermo`. So the primary deck could not resolve
electron thermo regardless of `u0`/`u1`: a deck-config gap stacked on top of the actual Wall A. The
prompt's Wall A ("the u1 electron misses the u0 DB entry") is only exercisable with an electron
thermo library loaded, which the **secondary** deck (`input_secondary.py`) does. To drive the
prompt's actual Wall A, `electrocatThermo` was added to `input.py` (keeping the `u1` electron); the
change is commented inline. With the electron thermo library present and the `u1` electron, both
fixes let the driver run through Wall A and Wall B.

## Where the pipeline stops now — Wall D (a successful outcome)

After the fixes, `python rmg.py docs/m7-preflight/input.py` runs through thermo (electron resolved,
renamed `[Lip]`), past the `Trange` access, into model generation — where the neutral CH₃Li generates
the lithium cation via reverse `Cation_R_Recombination`:

```
Created 1 new edge reactions
    [Lip](3) + [CH3](4) <=> CH3Li(2)
```

and then stops at **Wall D**, the Chemkin export's element-`E` balance guard, exactly as the ticket
predicts:

```
File "rmgpy/chemkin.pyx", line 1828, in rmgpy.chemkin.write_reaction_string
File "rmgpy/electron_balance.py", line 265, in check_electron_balance
  raise MechanismWriterError(
rmgpy.exceptions.MechanismWriterError: Reaction [Lip](3) + [CH3](4) <=> CH3Li(2) does not balance
  in the E pseudo-element and would be exported wrong: the equation
  "[Lip](3)+[CH3](4)<=>CH3Li(2)+e-(1)" has E=-1 on the left and E=1 on the right (reaction.electrons=1).
```

This is the loud, correct refusal the ticket says is a successful result — the E-balance guard
refusing a reaction the representation cannot express. It fired via the driver (not an isolated
probe), which means the pipeline now runs far enough to be honest about where it stops. Wall C
(reactor-side electron placement) was not reached because the Chemkin export in the initial
enlarge/save step raises Wall D first. Both are out of scope and untouched (`git diff` shows
`chemkin.pyx` and `electron_placement.py` unmodified).
