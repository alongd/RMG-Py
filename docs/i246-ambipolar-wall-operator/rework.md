# I-246 rework: the wall operator no longer refuses the chemistry it serves

This is the record of the rework asked for on top of the original I-246 wall operator (whose report
is `report.md`). It replaces the formula-keyed recycling refusal with a physically-grounded
ground-state rule, closes three validation gaps, states the regime of validity where a user sees it,
and enforces the single-cation support the internal report only admitted in prose.

All evidence is in `evidence/`. Every claim below is reproducible with the harness there against the
built module; the `.so` was proven by value after each rebuild (`strings` for a token added to the
`.pyx` -- `no usable thermochemistry`, `does not support anions`, `degeneracy threshold`).

## The headline: ion -> neutral resolution rests on the ground state, not the formula

**Decision.** An ion neutralised at a wall returns as the **electronic ground state** of its heavy
composition. Physically, wall neutralisation of a positive ion is a surface Auger process: the ion
captures a wall electron and the ~15.8 eV of ionisation energy is carried off by an ejected
(Auger) electron or into the surface, leaving the **ground-state** neutral. This is not the
electron-impact excitation that populates the Ar(4s) metastables at 11.5-11.7 eV; those are made in
the volume by electron collisions, not at the wall. So a wall that returns Ar+ as Ar* would be
inventing an excitation the surface process does not supply.

- Hagstrum, *Phys. Rev.* **96**, 336 (1954), "Auger ejection of electrons from tungsten by noble
  gas ions" -- Auger neutralisation of noble-gas ions at a surface yields the ground-state atom.
- Lieberman & Lichtenberg, *Principles of Plasma Discharges and Materials Processing*, 2nd ed.
  (Wiley, 2005), ch. 6 & 9 -- ambipolar wall loss and ion neutralisation at the boundary.

**Rule as implemented** (`_resolve_ground_state_neutral`, `rmgpy/solver/plasma.pyx`). Among the
neutral core species sharing the ion's heavy composition, the ground state is the one of **lowest
formation enthalpy** at the gas temperature -- which is the *definition* of the ground state and the
only discriminator that generalises. Multiplicity does **not** generalise: O2's ground state is the
triplet and its lowest metastable the singlet, so a "lowest multiplicity" rule returns the excited
state. That trap is recorded in the code comment so the next reader does not re-propose it.

Three cases the rule has to answer, and how it does:

| case | outcome |
|------|---------|
| two neutral states, one ground one excited (Ar, Ar*) | lowest enthalpy -> **ground**, logged |
| a molecular ion whose neutral has isomers, or any near-degenerate pair | **refuse**, name the declaration |
| excited-only deck (no ground present) | **not inferable** -- see the named gap below; the declaration resolves it |

**The precondition that had to be a refusal, not a fall-through.** The energy comparison needs
thermochemistry. Production core species carry it; a bare fixture (and, in principle, any deck built
without thermo) does not. The trap the ticket flagged is that a missing-thermo comparison must not
silently drop back onto the `len(matches)==1` "single neutral, use it" branch -- that is the
silent-wrong-answer path wearing a new hat. It refuses instead, naming thermo
(`evidence/before.log` H2: old message says only "ambiguous"; `after.log` H2:
`no usable thermochemistry`).

**The near-degeneracy threshold is a number with a reason: `k_B * T_gas`.** Two same-composition
neutrals closer in formation enthalpy than the thermal energy `k_B*T_gas` are appreciably
co-populated at the gas temperature -- the Boltzmann factor `exp(-DeltaE/k_B T_gas)` exceeds `1/e` --
so neither is "the ground state" in any operational sense, and picking one would move population
silently. It is temperature-aware and uses a quantity already in the reactor. Demonstrated
(`evidence/after.log`): a pair at `2*k_B*T_gas` **resolves** to the ground state (H3); a pair at
`0.5*k_B*T_gas` is **refused** naming the threshold and the declaration (H4). At 298 K the threshold
is 2479 J/mol; the ground/metastable argon gap is 1.11e6 J/mol, ~450x above it, so the deliverable
deck never comes near the boundary.

**The declaration.** `wallNeutralizationProducts={'Ar+': 'Ar'}` (`rmgpy/rmg/input.py`, and
`wall_neutralization_products` on the reactor) is the named escape for what energy cannot infer. It
is validated for shape and label existence at parse time and for charge/composition at reactor init;
a declaration naming a neutral absent from the core is refused by name (`after.log` H5, and the
`missing_neutral` test).

**Demonstrated on the actual collision** (`evidence/conservation.log`, and
`test_wall_resolves_metastable_deck_to_ground_state`): core = ground Ar, metastable Ar*, Ar+, e-,
`wall_recycling=1.0`. It initialises (it refused before), Ar+ recycles to **ground Ar**, and:

- **charge conservation:** net charge rate across the wall = 0 (exact) -- the common `nu` removes
  charge at rate `nu*(net charge)`, zero in a neutral gas;
- **particle conservation:** at gamma=1 the net heavy rate = 0 (every Ar+ returns as ground Ar, the
  metastable takes no wall flux); at gamma=0.5 the heavy loss is exactly `-(1-gamma)*nu*y_Ar+`.

## The three that had to close

**1. Neutral floor no longer silently accepted.** `compute_nu_wall` clamps neutral moles to
`wall_neutral_floor` to keep a Newton trial finite; `check_wall_support` used to validate only
positivity and the ionisation-degree ratio, so a state whose inventory had collapsed to the floor
(low `n_e/n_neutral`, tiny `n_neutral`) was **accepted** while `nu_wall` came off the clamp. The
ceiling bounds the *ratio*, not the absolute inventory -- so the old comment claiming the floor
"can never bind on a state the domain check would accept" was false, and is corrected in the source.
`check_wall_support` now refuses a neutral inventory at or below the floor (`evidence/before.log` F:
accepted; `after.log` F: `numerical floor`). The `test_check_wall_support_refuses_subfloor_inventory`
test pins it.

**2. Regime of validity is now in the user docs.** `documentation/source/users/rmg/input.rst` gains
a `.. warning::` next to the wall parameters stating plainly that this is one bulk ambipolar
eigenvalue with zero net current and no sheath -- defensible for an electropositive, unmagnetised,
diffusion-limited discharge (ICP or DC positive column, a few torr) and **not** reliable for a
high-pressure capacitive discharge, an electronegative gas, or a magnetised plasma.

**3. Test database path resolved from settings.** `test/rmgpy/solver/plasmaWallTest.py` no longer
hard-codes `/home/alon/Code/RMG-database-plasma/...`; it resolves `voronov.yaml` from
`settings['database.directory']` the way the rest of the suite does.

## Lower-priority, done

- **Single-cation support enforced** (`_resolve_wall_state`): a core carrying an anion, a second
  cation species, or a multiply-charged ion is refused with a named error, because the single
  `ionReducedMobility` cannot represent their transport and an anion is confined by the sheath, not
  lost at the wall (`before.log`/`after.log` A). Tests: `test_wall_refuses_anion_in_core`.
- **Non-finite electron population refused** in `check_wall_support` (`before.log`/`after.log` E) --
  `nan > ceiling` is False, so NaN/negative electron counts used to pass.
- **`diffusionLength` dimension checked** in `rmgpy/rmg/input.py`: `(2, 's')` was silently taken as
  `(2, 'm')` on the direct-length path (the shape path already checked units). Now refused; the
  valid `(2.03, 'cm')` spelling still resolves. Test:
  `test_direct_diffusion_length_dimension_is_checked`.
- **The "cannot fail" bitwise test** is rewritten
  (`test_wall_operator_residual_inert_and_jacobian_is_the_wall_linearization`). The old one compared
  two *wall-less* reactors, so it never evaluated the wall path -- vacuous. The rewrite pins two real
  invariants at an all-neutral state: the **residual** is bit-for-bit inert (the wall adds exactly
  zero), and the **Jacobian** difference from a wall-less reactor equals the analytic wall
  linearization (`-nu` on each charged diagonal, `+gamma*nu` into each recycle target). Note the
  Jacobian is *not* inert -- `d(nu*y)/dy = nu` even where `y=0` -- and claiming it were would be
  wrong; the invariant pins the value instead. Non-vacuity is proved in
  `evidence/bitwise_mutation.log`: a one-line corruption of `_apply_wall_terms` (`+1e-3` on the
  loss) makes the rewrite **fail** while the old two-wall-less pattern stays **green**.

## What I could NOT reach (named gaps)

1. **Excited-only deck is not auto-refused.** A core carrying *only* a metastable neutral (no ground
   state) presents Ar+ with a single composition-match, indistinguishable from an ordinary
   single-neutral deck. No thermo comparison can certify a lone metastable as "excited" without a
   ground reference present, and RMG carries no intrinsic ground/excited marker. So this deck falls
   on the `len(matches)==1` branch and uses the metastable silently. The ticket asked for a clear
   refusal here; the honest answer is that it is **non-inferable**, which is why the design carries a
   declaration. `wallNeutralizationProducts` is the mechanism -- and a declaration pointing at a
   ground state absent from the core *is* refused by name -- but absent any declaration the bare
   excited-only deck is used, not refused. A realistic metastable deck also carries the ground state
   (you do not add Ar* without Ar), so this is a corner, but it is a real one and I did not close it.

2. **Direct `PlasmaReactor(...)` construction still skips dimension checks** on `diffusion_length`,
   `ion_reduced_mobility`, `mobility_reference_density` and `ionisation_source` -- `_configure_wall`
   validates finiteness and sign but not dimension. The user-facing path (`plasmaReactor(...)` ->
   `_plasma_wall_kwargs`) checks all of them, and I added the missing `diffusionLength` check there,
   so a deck cannot hit this; only code constructing the reactor object directly can. I chose to
   report it rather than push `quantities`-based dimension machinery into `plasma.pyx` for an
   internal API, to keep the change inside the transport code minimal.

3. **`check_wall_support` does not assert full quasineutrality.** It now rejects a non-finite or
   negative electron population, but a finite, badly non-neutral state (many ions, few electrons)
   still passes the ionisation-degree ceiling. The wall's common-`nu` design conserves net charge
   regardless, and initial-composition neutrality is handled separately by the reactor, so I did not
   add a net-charge assertion here; flagged so it is a decision, not an oversight.

## Files touched

`rmgpy/solver/plasma.pyx`, `rmgpy/rmg/input.py`,
`documentation/source/users/rmg/input.rst`, `test/rmgpy/solver/plasmaWallTest.py`, and this
`docs/i246-ambipolar-wall-operator/` directory. No file under `rmgpy/molecule/`, `rmgpy/kinetics/`,
`rmgpy/data/` or the database was touched. Nothing pushed, merged or rebased.
