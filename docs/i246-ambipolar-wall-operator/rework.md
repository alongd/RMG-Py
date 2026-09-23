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

## Round 79: the class attacked, five HIGH and two MEDIUM closed

An adversarial pass found that the earlier rework had *narrowed* formula-matching, not replaced it.
The recurring root is one defect **class**: a species identified by a key coarser than the physics
distinguishes. It had produced three distinct defects -- ion→neutral recycle, neutral→cation for the
source, and (on a sibling branch) a transport estimator saturating metastable argon into ArH2. Rather
than patch a third site, the whole `plasma.pyx` file was audited for composition→species maps.

**The correct key is the heavy-atom skeleton: the standard InChI truncated before its charge (`/q`,
`/p`) layers** (`_skeleton_key`). Proven to sit exactly on the physics boundary: it unifies electronic
states (Ar and Ar* both `InChI=1S/Ar`, so the ground-state enthalpy tiebreak still applies) yet
separates constitutional isomers (DME `.../c1-3-2/...` vs ethanol `.../c1-2-3/...`, so a molecular ion
cannot be transmuted into an isomeric neutral). Element count is too coarse; `is_isomorphic` and the
augmented InChI are too fine. Both the ion→neutral and neutral→cation maps now key on the skeleton
(`evidence/round79_demo.log` §1-2: DME+ recycles to DME, not the lower-enthalpy ethanol).

Where even the skeleton is ambiguous the correct key is **user-supplied identity** via
`wallNeutralizationProducts` -- a design conclusion, not a gap. The excited-only deck (two argon
states, no ground) remains non-inferable and is reaffirmed, not newly closed.

Closed, each reproduced RED first (`evidence/round79_before.log`) then GREEN on the rebuilt module
(`round79_after.log`), with durable tests:

- **HIGH 1** -- isomer transmutation (skeleton key); and a non-finite (NaN/inf) enthalpy is now
  *unusable*, routed to the same refusal as absent thermo, not sorted silently into place.
- **HIGH 2** -- "exactly one cation" now rejects **zero** as well as two (an electron with no cation
  is not floating-wall physics); and the ionisation-degree ceiling gates on the **charged inventory**
  (`max(n_e, Σ positive-ion moles)`), not n_e alone, so 1% ions with 1 ppm electrons no longer slips
  under it.
- **HIGH 3** -- the external source is apportioned over **ionisable** neutrals only, so a non-ionisable
  bath gas (He) no longer sits in the denominator and swallows half the declared source
  (`round79_demo.log` §3: delivered/declared ratio 0.5 → 1.0000). The Jacobian mirror was changed in
  lockstep and checked against finite differences on the mixture path.
- **HIGH 4** -- the input writer now emits `wallNeutralizationProducts`, so the declaration survives a
  save/reload. The wall serialisation was extracted to `_format_plasma_wall` to make the round trip
  unit-testable.
- **HIGH 5** -- the wall-flux/wall-energy interface now exists: latched at accepted states only, with a
  machine-readable availability map. Full design in `interface-design.md`.
- **MEDIUM** -- refusal messages now name the input-deck keyword `wallNeutralizationProducts` a user can
  type, not the Python attribute; direct `PlasmaReactor(...)` construction now checks the *dimensions*
  of `diffusion_length` and `ion_reduced_mobility`, closing the last of named gap #2.

**Sibling ArH2 ticket (fourth instance of the class, worked separately).** A transport estimator
saturates a u2 (metastable) species into a molecule no atom type owns, because a family-level
containment cannot see a non-family consumer -- the same projection error, a different projection
(family membership rather than element count). The rule here applies there unchanged: key at the
finest identity the physics distinguishes. The skeleton key (charge/electronic-state-independent
InChI) is directly reusable as that identity if the two branches are reconciled; they should not
diverge on what "the same species" means.

## Round 83: the guard class, five HIGH and three MEDIUM closed

A second adversarial pass named the class from the other side. Round 79's class was the *identity
key*; round 83's is the **guard**: *a guard that reads a nearby quantity instead of the governing
one*. Before fixing the findings individually, every guard in `plasma.pyx` was audited for that
gap. The table is the durable deliverable -- it outlives the five repairs it organises:

| site | reads | physics gates on | verdict |
|------|-------|------------------|---------|
| `step()` diagnostics | *(latched nothing; only `advance` did)* | published diag time ≥ solver time on the **live** `simulate→step` path | HIGH 1 -- wrong path |
| `_skeleton_key` identity | standard InChI (mobile-H, tautomer-**merged**) | constitutional identity incl. H connectivity | HIGH 2 -- projection coarser than constitution |
| ion→neutral resolution | `len(matches)`; lowest enthalpy | true electronic **ground** (absolute anchor) | HIGH 3 -- relative order ≠ ground |
| no-cation guard | cation species **presence** | cation **inventory** vs electron inventory at the state | HIGH 4 -- presence not inventory |
| source deposit | silently skips when `y_ionisable≤0` | supported-neutral **inventory** must exist to receive the source | HIGH 5 -- silent zero, no accepted-state refusal |
| `mobility_reference_density` | magnitude only | number-density **dimension** (m⁻³) | MED -- value not dimension |
| `ionisation_source` | magnitude only | rate-density **dimension** (m⁻³ s⁻¹) | MED -- value not dimension |
| wall params when `has_wall=False` | *(nothing)* | wall-only options require a wall | MED -- accepted then ignored |
| `wall_neutralization_products` key | key/val are **strings** | key is an actual **cation species** | MED -- typo falls back to inference |
| `check_wall_support` neutral | **aggregate** `Σ y_neutral` | **per-species** non-negativity | MED -- aggregate not per-species |

One sentence covers column four: **presence where inventory governs (HIGH 1, HIGH 4); aggregate
where per-species governs (MED); a projection where constitution governs (HIGH 2); relative where
absolute governs (HIGH 3); residual-scratch on the `advance` path where the live `step` path
governs (the latch).** The audit found no sixth divergence: the `source_cation` ambiguity, the
orphan-ion check, and the alpha ceiling all read the governing quantity already.

**HIGH 2 and HIGH 3 collapse to one fact and one fix.** When 2+ neutral candidates share an ion's
skeleton, nothing in the deck certifies which is the ground-state product: they may be electronic
states (Ar/Ar*, where the lowest is the ground *only if the ground is present* -- an excited-only
deck is indistinguishable) or tautomers standard InChI merges (2-pyridone / 2-hydroxypyridine,
`evidence/round83_key_probe.log`), which are not electronic states at all. The FixedH InChI layer
separates tautomers on neutrals but does **not** compose with the charge-truncation the key needs:
on a cation InChI puts `/p`,`/f` before `/q` and shifts the formula, so truncating for
charge-independence strips the `/f` that distinguishes tautomers. A fourth projection was declined.
The fix (chosen with the reviewer) is **multiplicity requires a declaration**: 2+ candidates and no
`wallNeutralizationProducts` → refuse; do not pick. `_resolve_ground_state_neutral` -- the
lowest-enthalpy / k_B·T degeneracy picker -- was **removed**. Enthalpy leaves identity entirely and
survives only in the wall-energy interface, as a measurement of a quantity rather than a vote about
an identity. Cost, accepted as positively correct: a deck carrying two argon electronic states now
declares `wallNeutralizationProducts={'Ar+': 'Ar'}` -- one line in a deck that already declares
seven wall parameters, and a modelling claim the modeller should state.

The **single-candidate floor is named, not hidden** (docstring, `input.rst`): one neutral per
skeleton is returned because there is no alternative, not because it is certified as the ground
state; if the only neutral present is a different tautomer or an excited state, the ion returns as
that. No rule can distinguish that from one candidate; the declaration overrides it.

Closed, each reproduced RED first (`evidence/round83_before.log`, 14 red) then GREEN on the rebuilt
module (`round83_after.log`, 60 wall / 57 plasma / 101 input):

- **HIGH 1** -- `step()` (the entry `ReactionSystem.simulate` actually calls) now latches the wall
  interface after `check_wall_support`, so production runs publish live diagnostics, not the t=0
  latch. Test drives a real `step()` and asserts the published diagnostic time tracks the solver
  time -- the check that would have caught it.
- **HIGH 2 + HIGH 3** -- multiplicity requires a declaration (above); the energy picker removed.
- **HIGH 4** -- `check_wall_support` refuses a net charge outside quasineutrality (`|n_ion − n_e|`
  bounded by the existing net-charge tolerance), so a declared cation present at **zero moles** with
  electrons is caught -- the inventory the topology guard cannot see. Ordered after the alpha
  ceiling, so an over-ceiling non-neutral state still gets the more specific ceiling diagnosis. This
  closes deferred gap #3.
- **HIGH 5** -- `check_wall_support` refuses a declared source with no ionisable-neutral inventory
  left to receive it (only unsupported bath gas remaining), rather than depositing zero silently.
- **MEDIUM** -- direct-construction dimension checks on `mobility_reference_density` (m⁻³) and
  `ionisation_source` (m⁻³ s⁻¹); wall-only options supplied without a wall are refused rather than
  stored and ignored (`__reduce__` passes the wall-only densities back as `None` when wall-less, so
  a wall-less reactor still round-trips); and a `wallNeutralizationProducts` key naming no cation is
  refused rather than silently ignored; and per-species (not aggregate) neutral non-negativity.

## Round 88: the guard class again, four HIGH and two MEDIUM closed

A third pass, same class as round 83 for three of the four HIGH -- *a guard, a key, or a diagnostic
reading a quantity adjacent to the one the physics governs* -- plus one that is the round-83 latch
defect's structural twin: a term that moves the composition but is invisible to the code that asks
whether the composition can still move.

| site | reads | physics gates on | verdict |
|------|-------|------------------|---------|
| `_skeleton_key` identity | InChI **truncated at** the first `/q`,`/p` | nuclei + connectivity, one rule for **all charge states** | HIGH 1 -- `/i` and stereo sit AFTER `/q`, so truncation keyed a charged species by a different rule than a neutral: a ¹³C cation keyed as ¹²C neutral (transmutation) |
| `simulate` inert / termination | `char_rate` (gas-phase **chemistry** diagnostic) | the reactor's **total** flux | HIGH 2 -- a wall-only deck carries no `char_rate` but the wall moves it; declared inert at t=0 |
| quasineutrality bound (3 sites) | net charge vs **absolute** `1e-12 mol` floor | the **relative** imbalance `|net|/magnitude` | HIGH 3 -- a wholly-unpaired but small electron inventory (1e-13 mol) passes an absolute mole floor |
| `wall_neutralization_energy_flux` | `Σ ΔH · γ · (loss)`, gated on `γ>0` | energy owed per **ion lost** (γ governs mass, not energy) | HIGH 4 -- a pumping wall (γ=0) reports 0 W and labels it `available` |
| `_resolve_declared_neutral` | **first** label match in `matches` | the declared label must identify **one** species | MED -- a duplicate neutral label lets core ordering pick the product |
| `quasineutral_electron` flag | `bool(value)` (truthiness) | the flag's **value** | MED -- the string `'False'` is truthy, so it enabled the mode |

**HIGH 1 -- the key was two rules.** InChI orders its layers formula / c / h / `q` / `p` / stereo /
`i` (isotope) / …, so the isotope and stereo layers sit *after* the charge layers. Truncating at the
first `/q` or `/p` dropped them from a **charged** species and kept them on a **neutral** -- the one
asymmetry an ion→neutral map cannot have. Measured (`evidence/round88_high1_key_probe.log`, and
`test_skeleton_key_keeps_isotopes_...`): `[13CH3][O+][CH3]` → `.../h1-2H3/q+1/i1+1`, truncated to
`.../h1-2H3` = ordinary-carbon `COC`; the wall recycled ¹³C into ¹²C and the ionisation source did the
reverse -- nuclei transmuted and fabricated. And the escape hatch failed with it: the isotopic neutral
kept its `/i` (no `/q` to cut at), so it never entered the ion's candidate `matches` list and a correct
declaration naming it was *refused*. Fix, argued and adopted: wall neutralisation is a charge transfer
that conserves nuclei and connectivity but not charge or electronic state, so the key **removes the
`/q` and `/p` layers and keeps every other layer** (a first-char filter on the `/`-split tokens),
rather than truncating at them. Ar and Ar* still coincide (no layer between them); ¹³C and ¹²C now
separate on both charge states; the declaration resolves. This is one rule for every species
regardless of charge -- the actual defect.

**HIGH 2 -- the wall is invisible to "can this change?".** `res` carries the wall and source terms,
but `core_species_rates` deliberately does **not** (it is the gas-phase chemistry diagnostic that
model enlargement and the rate-ratio criteria compare edge fluxes against -- the comment guarding
that is correct and was kept). `base.pyx` then reads `char_rate ≤ floor` as "the composition cannot
change". On a wall-only deck (`simulate()`, gamma=0, no chemistry) the wall removes ~1.9e-2 mol/s of
ion-electron pairs while `char_rate = 0`, so the run terminated at t=0 with a warning that the
composition cannot change -- false. Two correctness criteria, each right, jointly impossible: "a
transport term is not a reaction flux" and "`char_rate` is the total flux". Fix: a polymorphic seam,
`ReactionSystem.get_non_chemical_char_rate()` returning `0.0`, overridden in `PlasmaReactor` to
evaluate exactly the terms `_apply_wall_terms` adds (the single source of the wall arithmetic, never a
re-derivation) at the accepted state, divided by V into `core_species_rates` units. The two inert
tests gate on `total_char_rate = √(char_rate² + non_chemical²)`; `char_rate` itself is untouched, so
the enlargement ratios and logs are unchanged, and a reactor with no non-chemical terms behaves
exactly as before (`total_char_rate == char_rate`). Demonstrated on `simulate()`
(`evidence/round88_high2_simulate.log`): the wall-only deck now integrates, the plasma decays from
x_ion=1e-4 to ~2e-17 over 0.16 s, and terminates at a genuine steady state; a wall-less inert deck
still stops at t=0.

**HIGH 3 -- quasineutrality is a ratio.** All three net-charge checks (`check_wall_support`,
`set_initial_conditions` under the algebraic electron, and the initial-composition warning) compared
`|net|` to `max(1e-12 mol, 1e-6·magnitude)`, which admits a state below *either* bound. On a small
charged inventory the 1e-12 mol floor was the larger term and swallowed a 100%-imbalanced state: 1e-13
mol of electrons with **no ion partner** reads as net −1e-13 mol < 1e-12, so it initialised. An
absolute mole tolerance cannot express a ratio -- the same imbalance passes or fails depending only on
how many moles the deck carries. Fix: the pure relative bound `|net| > 1e-6·magnitude`, scale-invariant
(100% imbalance refused whether the inventory is 1e-13 or 1e6 mol; a genuinely neutral deck, `|net|`
at accumulation roundoff ≈ `N·2.22e-16·magnitude`, admitted at every scale; `magnitude=0` gives
`net=0`, admitted). The absolute floor was deleted.

**HIGH 4 -- energy is owed per ion, not per recycled atom.** `wall_recycling` (γ) is a *mass-return*
fraction. A pumped ion still reaches the wall, still recombines with an electron there, and still
deposits `H_ion − H_neutral` at the surface. The energy term was gated on `γ>0` and multiplied by γ,
so γ=0 reported 0 W and -- worse than the number -- labelled it `available`, defeating the round-79
availability map whose whole point is that an absent datum is machine-readably absent. Fix: γ scales
only the mass return; the neutralisation energy is `Σ ΔH·(ion loss)`, independent of γ. `ΔH` is NaN
when the product or its thermo is unknown, which now correctly marks the aggregate `unavailable`
rather than a confident zero. (LOW, folded in: the energy fields are marked `available` only after a
final finiteness check.)

Closed, each reproduced RED first (`evidence/round88_before.log`, 8 red) then GREEN on the rebuilt
module (`round88_after.log`; suites 69 wall / 57 plasma / 101 input, and 23 steady-state / 6
zero-flux / 6 base / 5 simple / 6 liquid / 9 surface unbroken by the `base.pyx` seam). Charge and
heavy-atom conservation re-measured with the metastable present at γ ∈ {1, 0.5, 0}
(`round88_conservation.log`): net wall current exactly 0, heavy-atom loss exactly `−(1−γ)·(ion loss)`.

**What I could NOT reach.** No molecular-tautomer stereo case was exercised end-to-end (RDKit refuses
InChI for many charged aromatics/amides, keying them `None` -- the same wall round 83 hit); the
isotope path is proven, the stereo path is proven only at the key level by construction, not on a
built stereocentre. `get_non_chemical_char_rate` reuses `_apply_wall_terms`, so it overwrites the
residual-scratch fields at the accepted state -- harmless (the next residual overwrites them before any
consumer reads them), but it means the seam is not side-effect-free; a future refactor that reads
`wall_loss_rates` between a step and the next residual would need to know that.

## Round 90

Two HIGH, two MEDIUM, a LOW. HIGH 1 is a feature that did not do the thing it exists for; the rest are
the round-83/88 defect classes again -- a guard reading a quantity adjacent to the one the physics
governs, and an availability flag less sceptical than the number it vouches for.

| site | read (defect) | governs (fix) |
|------|---------------|---------------|
| electron seed (`_validate_electron_state`, `set_initial_conditions`) | "n_e must be strictly positive" -- always | n_e may be **exactly zero when a zeroth-order `ionisation_source` is declared**: that source seeds the first electrons and the discharge ignites from neutral gas. Without a source the only production is n_e-proportional chemistry, so zero is a fixed point -- still refused. |
| neutral floor (`compute_nu_wall`, `check_wall_support`, `jacobian`) | neutral **moles** ≤ 1e-6·(initial moles) -- extensive, history-dependent | neutral **number density** ≤ a fixed fraction of the reference density -- intensive, inventory-independent |
| `wall_flux` availability (`_latch_wall_diagnostics`) | marked `available` unconditionally | `available` only when the flux is finite |
| `diffusion_length` (constructor) | finite and positive | its **square** must also be usable -- a sub-underflow Λ (Λ²→0) is refused by name, not left to raise a raw `ZeroDivisionError` mid-solve |
| `_coerce_bool_flag` | ended in `bool(value)` -- coerces any type by truthiness | refuses anything that is not a bool, `None`, or a recognised boolean string |
| rate ratios (`base.pyx`) | `core/edge/network_rates / char_rate` with `char_rate=0` on a wall-only run → 0/0 NaN | divide by 1.0 when `char_rate==0`; the ratios stay chemistry-relative (only the inert/termination gates read the total) |

**HIGH 1 -- the external source could not ignite neutral gas.** `ionisationSource` exists so a discharge
starts from a declared physical mechanism rather than a numerical seed, but two guards (the
initial-composition check and the packed-state check) demanded a strictly positive electron with no
exception for a declared source. Every source test in the suite carried a seed, so the advertised path
was untested *and* unreachable, and the docs claimed an ignition the code refused. **Decision (argued,
not silent): admit exactly zero electrons when a strictly positive `ionisation_source` is present;
keep the strict-positive requirement otherwise.** The discriminator is whether a zeroth-order electron
production term exists: the source rate does not depend on n_e, so `dn_e/dt > 0` at n_e = 0 and the
state leaves the origin; the gas-phase ionisation is ∝ n_e, so without a source n_e = 0 is a genuine
fixed point. Both guards were relaxed in lockstep; the acceptance is a zero-electron deck with a source
integrating through `simulate()` (`round90_after.log`, HIGH 1 GREEN), and `input.rst` now states the
rule in the same commit.

**HIGH 2 -- the neutral floor was extensive, gating an intensive law.** `nu_wall ∝ 1/n_neutral`, a number
**density**, but the floor was `1e-6·(initial neutral moles)`. Round 74 found a pumping wall could
deplete inventory below the floor; round 83 made the validator refuse that state -- but the refusal
inherited the floor's wrong dimension, so *acceptance depended on how the deck got there*. The reactor
is isobaric, so the neutral density is pinned near `P/kT` regardless of how many moles remain: the
reviewer's depleting run held `n_neutral ≈ 1.6e23 m⁻³` the whole way down while the *moles* fell to the
floor and were refused, and `nu_wall` differed 10× (18.5 vs 185) across that history boundary for one
intensive state. Fix: the floor is a **density**, a fixed fraction (`1e-10`) of the mobility reference
density -- intensive, inventory-independent, far below any density an ambipolar-diffusion discharge is
run at. Consequence carried through: pinning the density makes `nu_wall` a *constant* on the clamped
branch (no V dependence, unlike the old moles clamp), so the Jacobian drops **both** the `1/y_neutral`
and the `dV/V` terms there, not just the first -- verified by a full-matrix finite-difference at a
below-floor state. The physical collisionless-transition floor (mean free path ~ Λ) is M9/sheath scope,
named not implemented. Acceptance: a small-moles/normal-density state is accepted with `nu_wall`
unclamped, and the floor is provably independent of the initial inventory (two decks, identical floor).

**MEDIUM -- availability honesty and strict coercion.** (a) An extreme mobility makes `nu_wall`, hence
`wall_flux`, non-finite, yet it was marked `available`; the flag now reads `unavailable` unless the
flux is finite (the round-88 lesson applied to a second field). A sub-underflow diffusion length
(Λ=1e-200, Λ²→0) raised a raw `ZeroDivisionError` from the callback; it is now refused by name at
construction. (b) `_coerce_bool_flag` ended in `bool(value)`, so `2`, `0.5`, `NaN`, `object()` enabled
quasineutral mode and `[]`/`{}` disabled it; it now refuses anything not a bool, `None`, or a boolean
string.

**LOW -- a char_rate=0 divide.** `base.pyx` builds rate ratios by dividing by the *chemistry* `char_rate`;
on a supported wall-only run every core rate is exactly 0, so `char_rate=0` and the ratios were 0/0
NaN with an `invalid value encountered in divide` warning. The consumers were enumerated in the code:
the ratios (and the branching numbers) are chemistry-relative enlargement signals and keep reading
`char_rate`; only the inert/termination gates read the total (`total_char_rate`, round 88). The 0/0 is
guarded by dividing by 1.0 when `char_rate==0`, so no NaN is laundered through `argmax`.

Closed, each reproduced RED first (`evidence/round90_before.log`, 6 red on the built module) then GREEN
on the rebuilt module (`round90_after.log`). Suites: 78 wall / 57 plasma / 101 input green, and the
full reactor set (steady-state, zero-flux, base, simple, liquid, surface, reverse-reconstruction) --
301 passed, 1 skipped -- unbroken by the `base.pyx` change. Charge and heavy-atom conservation
re-measured with the metastable present at γ ∈ {1, 0.5, 0} (`round90_conservation.log`): net wall
current exactly 0, heavy-atom loss exactly `−(1−γ)·(ion loss)`, `nu_wall=185.17` cross-checking the
reviewer's 185.15.

**What I could NOT reach.** The neutral-density floor's *value* (`1e-10·n_ref`) is a numerical floor, not
a physical-validity edge: the rigorous floor is the collisionless transition (ion mean free path ≈ Λ),
which needs a momentum-transfer cross-section and is M9/sheath scope -- named, deliberately not
invented here. Because the reactor is isobaric, that density floor is effectively unreachable in normal
operation (density is pinned at `P/kT` unless the gas is driven over the α ceiling, which is separately
refused), so the `check_wall_support` density-floor refusal is now a guard against trial/hand-built
degeneracy rather than a state a production run reaches -- correct, but it means the refusal path is not
exercised by an end-to-end run, only by a constructed below-floor state. HIGH 1's zero-electron path is
proven on the integrated-electron formulation; the algebraic (`quasineutralElectron=True`) formulation
admits zero identically but was not driven end-to-end from zero.

## Round 93: ignition is a new dynamical regime — two HIGH, three MEDIUM, one LOW

Round 90 made ignition from zero electrons real. Round 93 is the consequence: the zero boundary is a
**new dynamical regime**, and mechanisms that were correct for a *seeded* run are wrong there. Every
premise below was probed on the built module before a line was written (`evidence/round93_probe.py`,
`round93_probe_timing.py`); the six repairs were reproduced red-first (`round93_before.log`) and are
green after (`round93_after.log`). The nu_wall cross-check holds: 185.14 at the operating point, 185.17
at α=1e-6, against the reviewer's 185.15.

| # | reads | governs |
|---|-------|---------|
| HIGH 1 | quasineutrality tested **relatively** at every scale | the algebraic charge row is enforced to an **absolute** accuracy; a relative test on a sub-`atol` inventory measures solver noise |
| HIGH 2 | the criterion arms on an empirical slope reaching `R≥1` | a **saturating-from-zero** trajectory has a slope bounded by 1, and its electron is sub-`atol`: the two the arm was built for are inverted |
| MED 1 | each wall input is finite on its own | their **combination** (`nu_wall`) can still be non-finite; a subnormal `Λ²` is a finite input and a lethal divisor |
| MED 2 | one reduced mobility against the **summed** neutral density | the mobility is defined for **one bath gas**; a mixture needs a per-gas (Blanc's-law) mobility the model does not carry |
| MED 3 | the density floor built from `mobility_reference_density` | the transport reads `mu0·Nref` as a **product**; the floor must be invariant under `(Nref·c, mu0/c)` |
| LOW | `sqrt(Σ res²)` as the non-chemical rate | a representable-but-small flux **underflows when squared** and reads as inert while a source is declared |

**HIGH 1 — algebraic-mode ignition, decided by the `atol` seam.** `quasineutralElectron=True` failed to
ignite: at `t=3e-15 s` the whole charged inventory is ~2e-33 mol (seventeen orders under `atol=1e-16`),
the algebraic charge row's absolute residual is ~2.6e-34, and the *relative* guard read that as an 11%
imbalance and refused. Measured trajectory (`round93_probe.py`): below `atol` the relative imbalance sits
at ~6%, and the instant the inventory clears `atol` it collapses to ≤2e-16 (machine epsilon) on the row's
own. So the guard stands down **only while `magnitude < atol`, and only in the algebraic mode** — where
the row is enforced absolutely. This is *not* the absolute floor round 88 removed: that floor was `1e-12`
mol and admitted a genuinely-unpaired `1e-13` mol electron in the *integrated* mode; both sit above
`atol`, so round 88's example stays refused and the integrated guard is byte-for-byte unchanged. The stand-
down reads the actual `atol_array`, not a constant, so a user's tolerance choice moves it.

**HIGH 2 — the campaign's goal, arriving and not being recognised.** At `S=1e5` the reactor reaches
`n_e = S/nu_wall = 540.12 m⁻³` exactly and holds it, yet reported `reached=False`. Two causes, both
measured. (a) The trajectory `n_e=(S/nu)(1-e^{-nu t})` has log-log slope `nu t/(e^{nu t}-1)`, **strictly
below 1** for all t>0 — the arm (`R≥1`) can never fire, because it was designed for *decaying* transients
whose slope is unbounded, and a saturating rise is the opposite shape. (b) The electron saturates at
3.3e-21 mol, far under `atol`, so the generic residual sees only the neutrals — which a weak discharge
never perturbs, giving `residual=0` from the first step (`round93_probe_timing.py`: `max residual = 0`,
`streak=65`, `armed=False`). Arming-only would then fire at ignition, not saturation. The fix supplies the
electron's **own** slope so firing waits for it to go flat (density = `S/nu_wall`, verified to 8 figures),
and arms on **`t·nu_wall ≥ 1`** — the identical `t/tau≥1` standard, evaluated from the relaxation time the
wall knows (`1/nu_wall`) rather than from the bounded empirical slope. Gated on a live discharge (source>0,
electrons present), so a model that never started still reports not-reached — verified in the same test.
Both hooks default to none/False in `base`, so every ordinary reactor is unchanged (all reactor suites:
308 pass / 1 skip).

**MED 1 — refuse, don't merely flag.** `mu0=1e308` (overflows `D_a`) and `Λ=1e-160` (`Λ²=1e-320`,
subnormal) both give `nu_wall=inf`. Round 90 marked the term *unavailable*; that is not refusing a state
that cannot be integrated. Now `nu_wall` at the reference density must be finite at construction, and `Λ²`
must be a **normal** double. Round 90's availability test still stands — it now latches a *hand-built*
non-finite state, since the extreme mobility is refused before it can be built.

**MED 2 — the option the reviewer offered that the design forbids.** "Enforce single bath gas" would
refuse any deck whose neutrals span more than one heavy skeleton — but that is *every* multi-species
plasma: an inert He diluent, an isomeric neutral the wall must not transmute into, multiple ionisable
co-reactants (the source-apportionment tests). Enforcement breaks the isomer-transmutation, diluent and
apportionment tests by construction. "Composition-dependent mobility" (Blanc's law) needs a reduced
mobility *per* bath gas, which the model does not carry. So neither offered option is viable; the honest
resolution is to **warn**, once, naming the gases and the approximation, and state it at the mobility
keyword in `input.rst`. The Ar/Ar* deliverable shares one skeleton and does not warn.

**MED 3 — round 90's HIGH 2 in a new coordinate.** The transport reads `mu0·Nref` as a product, so
`(Nref·c, mu0/c)` leaves every `nu_wall` bit-identical — but the floor was `FRACTION·Nref` and scaled by
`c`, so physically identical inputs were accepted or refused differently. The floor is now
`FRACTION·PLASMA_LOSCHMIDT`, a physical constant invariant under the reparameterisation (and unchanged for
the default deck, where `Nref` *is* Loschmidt). Verified invariant across `c ∈ {1, 1e3, 1e-3}`.

**LOW — a norm that loses a nonzero flux.** `get_non_chemical_char_rate` squared each term; a source
delivering `res ~ S/Na` per species underflows to exactly zero once `res < sqrt(DBL_MIN) ≈ 1.5e-154`, so a
declared-and-admitted source read as inert. Replaced with a max-scaled L2 norm, so a nonzero flux stays
nonzero at any scale.

### What I could NOT reach / chose not to do
- **MED 2 is a documented approximation, not a physics fix.** A genuine Ar/He mixture still runs on the
  Ar⁺-in-Ar mobility against the summed density; only a warning and the docs mark it. A composition-weighted
  mobility is a real feature needing per-gas transport data and is out of this rework's scope.
- **HIGH 2's firing time depends on the residual's tolerance.** The electron channel fires at `nu·t ≈ 23`
  for `tol=1e-8` (`e^{-23} ≈ 1e-10`), so the reported density is `S/nu` to ~10 figures; a looser tolerance
  reports a slightly-less-saturated density. This is the criterion's own tolerance semantics, not new.
- The `t·nu_wall ≥ 1` arm uses `nu_wall` at the current step; for a discharge whose `nu_wall` drifts with
  composition this is the instantaneous relaxation time, which is the right local reading but not a global
  guarantee about a wildly non-stationary `nu_wall`.

## Round 96: the steady-state report was unsound in both directions — two HIGH, three MEDIUM, one LOW

Round 93's physics was right — ignition works and the saturated state lands on `n_e = S/ν_wall` to
seven figures — but the machinery that **reports** steadiness failed in *both* directions, and it
was round 93's own additions that broke it. Every finding was reproduced RED on the built module
first (`evidence/round96_before.log`: 6 failing) then GREEN (`round96_after.log`), one rebuild
(`round96_build.log`, `.so` proven by value: `available-single-bath-approximation`, `run-time worst
case`, `volumetric molar injection`, `not still rising`).

**The one sentence the criterion now answers, in code (`termination.py`):** *the steady-state
criterion measures the COMPOSITION — the mole fraction of every resolved state variable and its
log-log slope `R = |d ln x / d ln t|` — and every residual folded into it, generic (neutrals) and
external (the sub-floor electron), is that same intensive slope; arming and flatness are evaluated
per channel on it.* Three of this branch's HIGHs have now been extensive/intensive mismatches (R90
floor on initial *moles*, R93 floor on *mobility parameterisation*, R96 external residual on
*electron moles*); stating the quantity and making every residual measure it is the through-line.

**HIGH 1 — a false positive (the failure arming existed to prevent).** `steady_state_external_armed`
set the *global* `armed` once `t·ν_wall ≥ 1`, so a slow neutral reaction that had not run through
its own timescale was declared stationary the instant the electron saturated. Probed
(`evidence/round96_discriminator.log`): a deck with a slow `Ar → X` drain (`k=1e-12`, X seeded at
1e-2, filling on a ~1e12 s timescale) reported `steady=True` at 9 s while X's generic residual was
`6.4e-10` and **rising** (`1e-11 → 2e-11 → … → 6.4e-10`, doubling with t). Arming is **per quantity**:
the electron's arm vouches only for the electron and may license the system only while the generic
residual is **not still rising** toward its own `R ≥ 1` arm. This is threshold-free — the
discriminator is the *sign of the trend*, not a magnitude floor between the two decks (which would
be exactly the tuned constant the reviewer warned against): the inert control's generic residual is
`0.0` (flat → licenses), the drifting neutral's is rising (→ blocks). Acceptance both ways in one
build: the drift deck reports NOT steady, the no-slow-chemistry control still reports steady.

**HIGH 2 — the converse, a valid stationary composition that never terminated.** The generic residual
is on mole *fractions*; `steady_state_external_residual` measured electron *moles*. On a pumped
(`γ=0`) discharge whose fractions go stationary while the absolute inventory shrinks, the moles
residual reads large motion where the fraction is flat and poisons the MAX fold, so the backstop
fired `steady=False` on a stationary state. Fixed by measuring `x_e = n_e / Σy`. Where the neutral
bath is fixed (isobaric) the two coincide, so round 93's saturating-from-zero recognition is
unchanged; only the drifting-inventory case is corrected.

**MEDIUM 1 — the guard evaluated the wrong expression.** The construction check computed `ν_wall` at
the reference density (where `mu_i = mu0`), but run time forms `mu_i = mu0·Nref/n_neutral`. `mu0=1e20`
with `Nref=1e308` leaves `mu0` finite while `mu0·Nref = inf`, so `ν_wall` is infinite for every real
state yet the reference-density proxy read finite and admitted it. Now evaluated at the exact run-time
expression's worst case (the neutral-density floor, where `mu_i` is largest).

**MEDIUM 2 — a warning is not an availability state.** The reviewer rejected round 93's warn-only. A
multi-skeleton neutral bath (Ar+He) is still not refused — refusing forbids every multi-species
plasma, the premise the round-93 probe inverted — but the single-bath approximation is now recorded
as a queryable **availability state**: `wall_energy_availability['wall_flux']` and the electron-energy
flux become `available-single-bath-approximation` rather than plain `available`, so a consumer reading
the latched fluxes sees the approximation without a log line. The downgrade touches only fields built
from `ν_wall`, and only where they were otherwise `available` — it never launders a genuinely
`unavailable` (NaN) field into a usable one. A single-skeleton bath (Ar, or Ar and its metastables)
stays `available`. The warning is kept for the human.

**MEDIUM 3 — persistence was a step count, not physical time.** The window counted three accepted
solver steps, so persistence depended on the integrator's step controller. Now the flat run must also
span at least one e-fold in time (`ln(t_now/t_flat_start) ≥ 1`), the native scale of a `d/d ln t`
criterion. A residual of *exactly* zero is exempt: a structurally frozen composition is unambiguously
steady, and requiring the extra span there would push a fully-pumped (`γ=0`) discharge — whose steady
state is `n_e → 0` — past the mole floor into solver-noise-negative territory that the wall's domain
guard then refuses. (That exemption is what keeps `test_wall_only_deck_integrates_past_t0` green.)

**LOW — a subnormal source that injects nothing.** A positive but subnormal `ionisation_source`
(`5e-324`, `1e-320`) — or any value whose volumetric molar rate `source/Na` underflows — read as
`source > 0` and switched off the zero-electron ignition guard while `source·V/Na` injected exactly
zero. Refused at construction; the threshold is on `source/Na` (the per-volume molar rate), so the
normal-but-underflowing sibling (`1e-290`) is caught too.

**Test-quality (both found true).** `test_reduce_enumerates_every_constructor_parameter` was a lexical
surrogate — it searched the whole `__reduce__` body, so a parameter named only in a comment or the
docstring would satisfy it; it now strips the docstring and comments and searches only the returned
reconstruction tuple. And the blanket "every red state here is banked" header over-claimed: it is now
scoped to the defect-reproduction tests, with the invariant/property checks (closed-form re-derivation,
conservation, scaling, parameterisation invariance) named as asserting a positive property directly
rather than reproducing a defect.

### What I could NOT reach / chose not to do (Round 96)

- **The HIGH-1/HIGH-2 discriminator is trend-based, and a trend on a noisy residual is a heuristic.**
  A neutral whose residual is monotonically rising is caught; one that has already peaked and is
  decaying through the sub-tolerance band is treated as settled (correctly). A pathological neutral
  whose residual oscillated across tolerance without a clean trend is not something I exercised. The
  probe showed the real decks are cleanly monotone (rising) or flat/zero, so the heuristic holds for
  them; I did not prove it for an adversarially-constructed non-monotone residual.
- **MEDIUM 2 is still a documented approximation, not a physics fix** — the availability state records
  the approximation rather than removing it; a composition-weighted (Blanc's-law) mobility needs
  per-gas transport data the model does not carry, and remains out of scope.
- **MEDIUM 3's e-fold span is a natural unit, not a tuned one, but it does shift termination times.**
  Ordinary reactors now terminate up to one e-fold later than before (verified: the full solver suite
  stays green). The zero-residual exemption is what prevents that shift from colliding with the
  electron-negativity guard on vanishing-species decks; a deck whose steady state is a species reaching
  a small-but-nonzero floor rather than exactly zero was not separately exercised.

## Round 100: the discriminator's mechanism was unsound though its behaviour was right — three HIGH, two MEDIUM, two LOW

Round 96's three behaviours were verified correct in all directions (slow drift → not steady, shrinking
inventory → steady, zero-seed ignition → steady, inert → not steady). Round 100 did not touch that
behaviour; it replaced the *machinery* underneath it, which was reaching the right verdicts by unsound
means. **The one quantity the criterion is about is unchanged — the composition's mole-fraction slope —
and every residual still measures it; what changed is how the criterion decides a channel has stopped
moving and how long "stopped" must hold.** Evidence: `evidence/round100_before.log` (9 red on HEAD),
`round100_after.log` (9 green on the rebuild), `round100_build.log` (clean `make clean && make build`,
`.so` proven by value: `wall_single_bath_approximation`, `run-time volumetric injection`,
`steady_state_relaxation_time`).

- **HIGH 1 — the external arm was a two-point comparison with a permanent latch.** It compared the
  generic residual to the single preceding sample and, once any non-increasing sample appeared (a noisy
  dip, two equal samples, or the very first step, whose predecessor is `nan`), armed the external
  channel *forever*. A still-rising neutral could arm the criterion. Fix: (a) "departing" is now a
  property of a **sequence** — the generic channel counts as still climbing until it has set no
  step-over-step rise for `window` consecutive samples, so one flat/noisy sample cannot license the arm;
  and (b) the external arm is **re-evaluated every step, not latched** — it is a conjunction with a live
  condition (no generic channel departing), so a neutral that settles and then resumes moving withdraws
  it. The generic `R >= 1` arm stays latched, because passing your own relaxation time is a historical
  fact, not a reversible condition. Threshold-free throughout: the discriminator is the sign of the
  trend over a window, never a magnitude floor.

- **HIGH 2 — persistence still depended on accepted solver steps, and an exact-zero waiver aliased
  oscillations.** `streak >= window` was mandatory, so fine stepping terminated inside a plateau that
  coarse stepping could not, and a residual of exactly zero waived the physical span entirely (equal
  endpoints do not prove a frozen structure). Fix: the binding requirement is a **physical span**
  (below); the step count survives only as a two-sample fluke guard (`streak >= 2`), deliberately **not
  sufficient**; and the **exact-zero waiver is removed** — a zero residual earns the same span
  confirmation as any other flat tail.

- **HIGH 3 — the physical window was anchored to absolute log time.** A tail first seen at `t0` had to
  survive to `e·t0`, i.e. another `1.718·t0` *regardless of the system's own relaxation time*, so the
  same physics converged or not depending on when the window happened to open; and the first flat
  interval was discarded (`_t_flat_start` was set to the interval's *end*). Fix: persistence is anchored
  to **one system relaxation time**, which the reactor supplies (`steady_state_relaxation_time` →
  `1/nu_wall` for a wall-bounded discharge, `nan` for an ordinary reactor, which then falls back to the
  absolute e-fold). `_t_flat_start` is now the flat interval's *first* endpoint (`t_prev`). Anchoring to
  τ made the confirmation short enough that removing HIGH 2's exact-zero waiver did **not** reintroduce
  the γ=0 collision (the absolute e-fold, when `t0` was large, forced far more than one τ of extra
  integration and drove `n_e` past the wall guard; one τ does not). So the electron-negativity guard was
  left untouched, as the contract requires — τ-anchoring, not a looser guard, is what resolves it.

- **MEDIUM 1 — the source guard checked a proxy.** `__init__` validated `source/Na`, but the residual
  injects `source * V / Na`. At an extreme-but-finite volume that product overflows or underflows while
  `source/Na` looks finite. The same shape as round 96's `nu_wall` overflow. Fix: the run-time
  expression is re-checked at the actual initial volume in `set_initial_conditions`.

- **MEDIUM 2 — the availability label documented the wrong transport rather than gating on it.** A
  multi-skeleton neutral bath was run silently on the single ion mobility, recorded only as a label a
  consumer might read. Round 96 rejected warn-only; round 100 rejects label-only. Fix: the mixture is
  **refused at construction unless the user opts in** with `wallSingleBathApproximation=True`
  (constructor arg + camelCase DSL keyword, round-tripped by `__reduce__` and `save_input_file`). Opt-in
  is the honest middle between refusing every multi-species plasma and running silently on transport that
  does not describe the gas; once opted in, the fluxes still carry `available-single-bath-approximation`.
  This is a judgement call the reviewer may still contest; the deferred fallback (a bath-gas-density
  `nu_wall` using only neutrals sharing the cation's skeleton) remains available but trades one
  approximation for another and was not done.

- **LOW — the tiny-source test asserted only a collapsed scalar**, which passes if only one member of
  the ion–electron pair is injected. Strengthened to assert the residual injects a positive rate into
  **both** the electron and its cation.

- **LOW / test-quality — the `_feed` harness stamped a clock onto `term._feed_step`**, an attribute the
  production object never sets, proving nothing about how the criterion is really called. The clock now
  belongs to the harness (keyed per term), and `term` is driven only through its real `update` API.

- **Test-quality — the red-state brief.** The claim that every test banks a red state is now scoped by
  an **observable marker**: a defect-reproduction names its round/finding in its docstring (and is
  reproduced red in that round's `before.log`); an invariant/property/acceptance check states a property
  and claims no red. A reader audits the docstring tag against the round log, not a separate manifest.

The through-line from rounds 90/93/96 (state one quantity, make every residual measure it) held: this
round did not add a quantity, it made the *decision procedure* over the existing quantity
step-controller-independent and threshold-free.

### The trend discriminator, again — found by probing, not reasoning

The HIGH-1 fix (like round 96's) rests on the sign of a trend over a window, not a magnitude floor. I
validated it by driving the built `TerminationSteadyState` directly (`high1_probe`): a neutral rising
throughout with one flat sample, electron past its relaxation time, must not arm; a settled neutral
must; a settled-then-departing neutral must un-arm. All four directions confirmed before the fix was
written. The permanent latch's starkest failure surfaced only under the probe — it armed on the *first*
step, because `nan`-predecessor reads as "not rising."

### What I could NOT reach (named gaps)

- **The residual is a scalar max over species.** HIGH 1's "not rising for `window` samples" reads the
  aggregate `r_gen`, so a slow neutral whose rise stays hidden beneath an earlier, larger species'
  decay is not separately detected. This is inherent to a scalar `d/dln t` criterion and unchanged from
  before; a per-species trend test is out of scope.
- **Coarse/fine and t0-independence are shown at the unit level**, driving `update()` with controlled
  schedules over the identical trajectory — more direct than hoping the integrator picks coarse vs fine
  steps, but not a demonstration on a real adaptive schedule. The integration tests exercise the real
  `simulate()` path for the four behaviours.
- **The two-sample guard can still cost coarse stepping one extra step** before it terminates (a single
  step spanning a whole τ has `streak == 1`); the verdict is unchanged, only the exact stop time shifts
  by one step. Reducing the guard below two would let a lone flat step terminate, which is the fluke it
  exists to reject.
- **MEDIUM 2 remains an opt-in to an approximation, not a physics fix** — the composition-weighted
  mobility needs per-gas transport data the model does not carry.

## Round 106: the discriminator got a slope but kept watching the aggregate — one HIGH closed, one HIGH rebutted, three MEDIUM, a census made executable

Round 100 replaced the external arm's two-point latch with a sequence-trend test — but the trend was
read off the **aggregate maximum** residual, the same scalar `r_gen` the round-100 "What I could NOT
reach" note flagged as out of scope. Round 106's reviewer made that gap a HIGH and was right to: a
guard reading a quantity adjacent to the one the physics governs is this campaign's oldest pattern,
and it had reappeared *inside* the fix for it. Evidence: `evidence/round106_before.log` (4 defect
reproductions red on HEAD, the census enforcement test already green), `round106_after.log` (all 5
green on the rebuild), `round106_build.log` (clean `make clean && make build`; `.so` proven by value:
`_source_total_at_volume`, `at volume V=... forms a run-time volumetric injection`).

- **HIGH 1 (closed) — the departing test read the aggregate maximum, not each species.** "Not rising
  for `window` samples" was evaluated on `r_gen = max_i |dln x_i/dln t|`. A species climbing toward its
  own transient is invisible whenever a *different*, decaying species holds the maximum — so the
  aggregate "stops rising", the system is judged no longer departing, and the external (electron) arm
  is let through while a live neutral is still moving. It also compared the slopes of two *different*
  species across any step in which the maximum changed hands. Fix: the slope analysis now returns a
  per-species dict (`_slope_analysis`, keyed by the core-species **integer index** — hashable and
  stable as the live set changes), and "departing" is tracked **per species** (`_slope_prev`,
  `_steps_since_rise`): each species counts as still climbing until it fails to rise against *its own*
  previous sample for `window` consecutive samples, and the generic channel is departing while ANY live
  species still is. The generic `R >= 1` arm stays on the aggregate maximum — that is exactly the right
  quantity for "did ANY species reach its relaxation time", and it is a latched historical fact.
  `compute_residual` is now a thin wrapper over `_slope_analysis`, preserving its `(r, label)` contract
  for every existing caller. Reproduced by `evidence/round106_high1_probe.py` (asserts + non-zero
  exit): the reviewer's hidden-rising sequence (`a`: 5e-8→1e-7→1.5e-7 rising; `b`: 6.77e-7→6.01e-7→
  5.26e-7 falling, always the max) does NOT arm, `worst_label` correctly names the falling `b`, and a
  genuinely-settled control still arms.

- **HIGH 2 (rebutted, with measurement) — the period-9 decade-aliasing case cannot arise.** The finding
  posits a periodic trajectory that arms at t=1 then reads exactly zero at t=10 and t=100, so an
  endpoints-only span check returns steady while the interior swings. Its load-bearing premise is that
  `update()` is sampled at **decade-spaced points**. Probed and inverted: `ReactionSystem.simulate`
  drives the integrator in pydas **intermediate one-step mode** (`self.step(step_time)`,
  `base.pyx:796`) and calls `update()` with `self.t` after **every internal DASSL step**.
  `evidence/round106_high2_cadence.py` (asserts + non-zero exit) instruments a genuinely relaxing wall
  reactor and measures the actual cadence: **415 samples, median 0.008 decade, max 0.477 decade, zero
  steps spanning a full decade.** DASSL's local-truncation-error control cannot accept a step that
  spans an oscillation period, so a flat residual over an accepted step is the integrator certifying
  the composition is smooth across it; a period-9 oscillation shows non-zero residual at the fine
  interior steps the criterion actually sees. An interior-variation guard would be redundant (the
  streak already requires every step's `R < tol`, which bounds the total variation across the span
  below `tol·Δln t`) and would fight DASSL's legitimate step growth in the flat tail — a false negative
  in the primary use case. Verifier item 6 invited this; no code guard was added, and the acceptance
  case is answered by showing the sampling it assumes does not occur.

- **MEDIUM (closed) — `window` was documented as the minimum accepted-sample span but hard-coded to
  two.** The termination test used `streak >= 2` regardless of the `window` a deck set, so window=3/10/
  100 all terminated after two flat samples — a knob that did nothing. Fix: `streak >= self.window`;
  the default and floor are now **2** (a flat interval needs two endpoints; a single step, however
  long, is not an interval), and raising `window` demands that many flat samples. It remains never
  *sufficient*: the physical span must also hold. Documented honestly in `input.rst` and reproduced by
  `test_the_window_parameter_is_honoured_as_the_minimum_sample_span` (window=4 fires on the 4th flat
  sample, not the 2nd). Four existing steady-state tests that had encoded the old hard-coded floor were
  moved to `window=2` so their two-sample narratives stay valid; the input default assertion moved from
  3 to 2.

- **MEDIUM (closed) — the source guard and the source computation evaluated different expressions, for
  the third time.** `set_initial_conditions` validated `source*V/Na` only at the **initial** volume,
  while the residual and Jacobian recompute `source*V/Na` at the **evolved and Newton-trial** volume. A
  source finite and positive at V0 overflows (or underflows to zero) once scaled by an extreme trial
  volume, and was then applied as an infinite (or vanishing) rate silently. Fix: a single validated
  helper `_source_total_at_volume(V)` (raising `PlasmaStateError` on a non-finite or subnormal product)
  is called from all three sites, so the guard travels with the value to every volume the solve reaches.
  Reproduced by `test_the_ionisation_source_is_validated_at_the_evolved_volume_not_only_the_initial_one`
  (a source valid at V0; a wild trial neutral amount inflates V until `source*V/Na` overflows; both
  `residual` and `jacobian` refuse it).

- **MEDIUM (closed) — the mixed-gas opt-in was truthiness-coerced.** `wall_single_bath_approximation`
  was set with `bool(value)`, so the string `"False"`, the int `2`, and `NaN` all opted in — and a deck
  writing `wallSingleBathApproximation="False"` silently enabled the approximation it meant to decline.
  The identical fix `quasineutralElectron` already carries: `_coerce_bool_flag` in the constructor
  (bool/None through, boolean-like strings by meaning, everything else refused by name), and `input.py`
  passes the raw value rather than `bool()`-casting it. Reproduced by
  `test_wall_single_bath_approximation_is_coerced_by_value_not_truthiness`.

- **Census (made executable) — the red-state brief was prose, and a hand-count disputed it.** The brief
  scopes evidence by an observable marker: a test is a banked defect reproduction iff its docstring
  names the round/finding it closes, otherwise it is an invariant/property check that asserts its
  property directly. That scoping is now enforced by
  `test_every_solver_test_is_a_tagged_defect_reproduction_or_asserts_a_property`, which walks the AST of
  both changed solver test files and fails, by name, any `test_*` that neither carries a finding tag nor
  contains an assertion (a `pytest.raises`/`warns` context counts). After the changes: **plasmaWallTest
  119→ of which 88 functions, 52 tagged, 87 asserting, 0 neither; steadyStateTest 31 functions, 11
  tagged, 31 asserting, 0 neither — 0 of 119 lack evidence.** The reviewer's 22/85 and 44/114 counted
  untagged invariant tests as gaps; the executable check confirms every one of them asserts a property
  (so it carries the evidence the brief requires) while catching any genuinely assertion-free untagged
  test in the future. The round-106 probes (`round106_high1_probe.py`, `round106_high2_cadence.py`)
  assert and `sys.exit(1)` on failure, unlike the print-only `round100_high1_probe.py`.

### What I could NOT reach (named gaps, round 106)

- **The per-species trend still reads a discrete two-sample slope**, so a species whose rise is slower
  than one sample interval registers as flat for that step; it is caught once its accumulated rise
  exceeds the sample-to-sample noise, not instantaneously. This is inherent to a finite-difference
  `d/dln t` and is why the arm waits for a *window* of non-rises, not one.
- **HIGH 2 is answered by rebuttal, not a guard.** If a future change moved `simulate` off intermediate
  one-step mode (e.g. stepping directly to decade output times), the cadence premise would flip and the
  aliasing case would become reachable; `round106_high2_cadence.py` is the tripwire that would catch
  that regression (it asserts the fine cadence and exits non-zero if a step ever spans a decade).
- **The source helper raises inside `residual`/`jacobian`.** A wild Newton trial that overflows
  `source*V/Na` now aborts the solve rather than letting DASSL reject the trial. This is the intended
  fail-loud behaviour (an infinite/vanishing source silently applied is the defect), and such a volume
  is astronomically unphysical, but it is a hard stop rather than a graceful trial rejection.

## Round 109: the poison is finite-checked away — one HIGH closed, three MEDIUM, a census of the whole class

The round-106 external-arm fold folded the reactor's electron residual into the generic residual by
MAX so a flat tail waits for the slowest of everything. But it wrote the MAX as
`if not np.isfinite(r) or external_residual > r: r = external_residual`, and `not np.isfinite(r)` is
true for **both** `nan` and `inf`. `_slope_analysis` emits `inf` when a species crosses UP through the
floor (appears from nothing) — "emphatically not steady" — so a small flat electron residual
**replaced** that poison instead of losing a MAX to it. The second half made it terminate rather than
merely mis-report: the per-species departing counters are advanced only on an `ok` step, so on an `inf`
step nothing was recounted or pruned, a composition settled a step earlier kept its counters,
`generic_departing` stayed `False`, `armed_external` stayed `True`, and the streak advanced on the
substituted residual. Under model enlargement a species crossing the floor is ordinary, so this was a
live path. No compiled code changed this round: the whole fix is in the pure-Python
`termination.py`. Evidence: `evidence/round109_before.stdout.log` (the two HIGH tests red on HEAD
`0cda42441`, the two MEDIUM property tests already green), `round109_after.stdout.log` (288 passed, 1
skipped on the fixed tree).

- **HIGH (closed) — a finite external residual erased the generic `inf`.** Two coordinated changes,
  both in `update()`. (1) The fold now reads `if np.isnan(r) or external_residual > r`: `np.isnan`
  singles out the *no-information* case a `nan` generic residual means, while an `inf` — *strong
  negative information* — is preserved (`external_residual > inf` is False), falls through to the
  non-finite streak guard, and correctly returns not-steady. (2) The counter update now branches on
  the analysis status explicitly, because `nan` and `inf` are **different answers**: an `ok` step
  counts each species' own trend and prunes; an `inf` step **discards** the per-species history
  (`_steps_since_rise.clear()`, `_slope_prev.clear()`) because the composition changed structurally and
  a stale "settled" count must not vouch for the new state; a `nan` step **leaves the counters exactly
  as they are** (no information observed). Both behaviours carry a code comment saying which is intended
  for which. Reproduced through the real criterion by moving a mole value across the real floor (never
  by injecting a non-finite): `test_a_species_appearing_while_the_external_channel_is_flat_does_not_return_steady`
  (the appearance returns not-steady and reports `inf`, while a byte-identical control with no
  appearance still terminates) and `test_nan_and_inf_steps_treat_the_departing_counters_differently`.

- **MEDIUM 1 (answered in docstring and user docs) — `window` gates two quantities, and they pull the
  SAME way.** The finding's premise, that raising `window` tightens persistence while *weakening* the
  departing test, is inverted. Raising `window` requires a longer flat streak (stricter persistence)
  AND holds each species as "departing" for more consecutive non-rises (later external arm) — both make
  the verdict strictly **more** conservative, so a larger `window` can only push termination later,
  never earlier. Splitting the knob would let a user buy a combination ("stricter streak, weaker
  departing test") that means nothing physical here. Documented in the `TerminationSteadyState`
  docstring and in `input.rst`, and pinned by
  `test_raising_window_makes_the_criterion_no_less_conservative_in_both_roles` (the firing step is
  monotonically non-decreasing across `window` = 2,3,4,5, driven through the external-arm gate so both
  roles are active).

- **MEDIUM 2 (answered by a test that grows the core) — index keys are safe because `reset()` wipes
  them between runs.** `_slope_prev`/`_steps_since_rise` are keyed by core-species integer index. Within
  one `simulate()` the core is fixed, so an index means the same species for the whole integration; the
  core only grows across `simulate()` boundaries (enlargement), and `simulate()` calls `reset()` first
  (`base.pyx:754`), clearing both dicts. `test_departing_counters_are_keyed_by_index_and_cleared_between_runs`
  demonstrates the hazard concretely — feeding an enlarged core WITHOUT a reset does inherit a stale
  count at a reused index — and then the guard: after `reset()` the dicts are empty, so no index can
  carry a meaning from a smaller core.

- **MEDIUM 3 (accept/reject enumerated, refusal reaches the user) — `_coerce_bool_flag`.**
  `test_coerce_bool_flag_accept_reject_set_is_enumerated_and_the_refusal_reaches_the_user` fixes the
  full contract: **accepted** — `True`/`False`/`None` (None → False), and the boolean-like strings
  `true/1/yes/on` → True, `false/0/no/off` and the **empty string** → False (case- and
  whitespace-insensitive); **rejected with a raised `PlasmaStateError`** — the motivating `"False"`'s
  sibling typos, any unrecognised string, `2`, `0`, `0.5`, `NaN`, `inf`, and any list/dict/object. The
  test also drives the rejection through the reactor constructor to show it propagates to the user, not
  into a log line.

### Census — every site where one channel's value meets another's, and what it does with `nan` vs `inf`

The defect is a **class**: a sentinel carrying "not steady / invalid" erased by a later stage's
finite/truthiness check. Enumerating every cross-channel combine / substitute / compare site in the
branch (a value from one channel — generic-neutral residual, external-electron residual, electron
population, cation population, differential Jacobian rows — meeting a value from another):

1. **`termination.py:260` — the fold** (generic residual vs external residual; substitute-or-MAX). The
   external side is gated by `np.isfinite(external_residual)`, so a `nan`/`inf` external is ignored.
   Generic side after this round: `nan` → substitute external (no information); `inf` → **preserved**
   (external cannot erase it); finite → MAX. **This was the one site with the erasure defect; now
   fixed.**
2. **`termination.py:334` — the external arm gate** `armed_external = external_armed and not
   generic_departing` (external vouch vs generic departing verdict). Both operands are booleans;
   `generic_departing` derives from integer counts and `external_armed` from the reactor. Non-finite
   slopes are dropped upstream in `_slope_analysis` (only finite slopes enter the dict), so no
   non-finite value can reach here. No erasure pathway.
3. **`termination.py:337` — the combined arm** `armed = armed_generic or armed_external` (two booleans).
   No non-finite pathway.
4. **`plasma.pyx:2704` — the quasineutrality residual row** `delta[electron_index] = charge_row_scale *
   _net_charge(y)` (electron population vs cation populations). A non-finite population makes the row
   non-finite, which the integrator's own error control rejects — `nan` and `inf` both propagate as a
   bad residual and are **not** erased into a spurious steady/valid state.
5. **`plasma.pyx:2731 (`_set_charge_row_scale`)** `scale = max(differential); charge_row_scale = scale
   if isfinite(scale) and scale > 0 else 1.0` (the charge row vs the differential rate rows). A
   non-finite `scale` (`nan` or `inf`) falls back to `1.0`. This is an **exact** row scaling that
   changes no solution, so the fallback cannot erase a steady-state sentinel; it only affects
   conditioning.
6. **`plasma.pyx:1141` — the initial net-charge diagnostic** (electron vs cations on the initial
   composition) degrades a non-finite `net` (`nan` or `inf`) to a warning and returns. This is an
   advisory neutrality check, not a steady-state or validity decision, and the electron population's
   own non-finiteness is separately **refused** at `plasma.pyx:2434` (`if not np.isfinite(e0) ...
   raise`). So the non-finite is caught by a refusal elsewhere; the diagnostic's graceful skip erases
   no sentinel.

**Count: 6 cross-channel sites. Sites where a sentinel meaning "not steady / invalid" could be erased by
a finite/truthiness check into a false positive: 1 before this round (the fold, site 1); 0 after.** The
external channel's own vocabulary is worth recording: `steady_state_external_residual`
(`plasma.pyx:2259`) returns **only** `nan` or a finite slope, never `inf` — so the fold's `inf` can only
ever arrive on the generic side, exactly where site 1 now handles it.

> **Round 110 correction.** The sentence above — "returns only `nan` or a finite slope, never `inf`" —
> was FALSE, and the round-110 HIGH 1 addendum proves it through the compiled hook: the hook guarded the
> electron MOLES (`ne_now > 0`) but took `log()` of the FRACTION `xe = ne/tot`, and a positive subnormal
> `ne` whose ratio underflows to `0.0` gave `log(0) = -inf` → `+inf`. It is fixed at the source (guard the
> fraction after the division; return `nan`), so the vocabulary claim is restored rather than merely
> asserted. The census count is also revised below: it omitted the FLUX/TRANSPORT channels and so was
> **6, not the exhaustive total — the true count is 8** (see the round-110 section).

This branch has now produced the same shape four times — `nu_wall` guard vs its computation, source at
initial vs evolved volume, aggregate vs per-species trend, and now a finite-check vs the specific
non-finite value that carries the meaning. The census above is the falsifiable form of "and nowhere
else".

### What I could NOT reach (named gaps, round 109)

- **A `nan` step still resets the flat streak** (via the shared non-finite guard) while leaving the
  per-species trend memory intact. This is deliberate — the two facts are independent (a no-information
  interval breaks *continuous* flatness but tells us nothing about any species' trend) — but it means a
  run peppered with degenerate intervals re-pays the streak each time even though no species moved. In
  practice DASSL does not emit repeated identical time points, so this is a latent property, not an
  observed cost.

## Files touched

Round 109 changed only `rmgpy/solver/termination.py`,
`documentation/source/users/rmg/input.rst`, `test/rmgpy/solver/steadyStateTest.py` and
`test/rmgpy/solver/plasmaWallTest.py` — no compiled code, so no rebuild. Across the whole rework:
`rmgpy/solver/plasma.pyx`, `rmgpy/solver/base.pyx`, `rmgpy/solver/base.pxd`,
`rmgpy/solver/termination.py`, `rmgpy/rmg/input.py`, `documentation/source/users/rmg/input.rst`,
`test/rmgpy/solver/plasmaWallTest.py`, `test/rmgpy/solver/steadyStateTest.py`,
`test/rmgpy/rmg/inputTest.py`, and this `docs/i246-ambipolar-wall-operator/` directory. No file under
`rmgpy/molecule/`, `rmgpy/kinetics/`, `rmgpy/data/` or the database was touched. Nothing pushed, merged
or rebased.

## Round 110: the value is ignored, the permission is not — two HIGH, one MEDIUM, a census whose detector was fixed

Closed against branch head `31487b60d` under the owner's binding ruling of 2026-09-23. The wall-transport
physics is **not** the subject and was not reopened: the measured wall-loss frequency is preserved
(re-measured below). Two HIGH, one MEDIUM, three test findings, and a census whose *detector* — not just
its list — was corrected.

### The numpy.bool_ type leak (owner's rebuild was RED here, and this is why)

The round-109-plus-partial branch, rebuilt, failed three `TerminationSteadyStateLatchTest` cases on
`assert term.armed_external is False` — an assertion that reads like `False is False`. The value was
correct; the **type** was not. The round-110 HIGH 1 arm fix added `np.isfinite(external_residual)` into
the `and`-chain, and `np.isfinite` returns a `numpy.bool_`; Python's short-circuiting `and` returns that
scalar when it is the deciding operand, so `armed_external` (and, through the `or`, `armed`) silently
became a `numpy.bool_` — truthy/falsy but failing an `is`-identity test a caller may make. This is the
same class as HIGH 1 itself: a value *almost* right on a contract that says *exactly* right. Fixed by
coercing both attributes with `bool(...)` at their one assignment site (`termination.py`). **What else
changed type the same way:** only `armed_external` (and `armed` transitively). `armed_generic` is assigned
a literal `True`, never a numpy expression; `residual` is a Python `float` (`external_residual` arrives as
a `cpdef double`, `r_gen` is `float(...)`-wrapped); `streak` is an `int`. The three failing tests'
*control* sub-cases (a settled channel that DOES arm) were fed `external_armed=True` with no
`external_residual`, so it defaulted to `nan`; under the round-110 contract a vouching electron channel
reports a finite slope (a sub-floor-but-positive electron gives a finite residual, not `nan`), so those
calls now pass `external_residual=1e-9`. That restores the exact generic-departing gating they test and
is faithful to the production path — not a relaxed assertion.

### HIGH 1 — a non-finite residual must never authorise anything

Round 109 made the fold *ignore* a non-finite `external_residual`; it left the boolean `external_armed`
(same call, independent) to grant the arm alone, and it dropped a non-finite external residual as
"criterion unavailable" while a flat generic channel went on to terminate. Reproduced through production
`update()`, `window=2`, flat generic channel: `external_residual` ∈ {`1e-9`, `nan`, `+inf`} all gave the
byte-identical `[False, True, True]`. Two-part fix, both in `termination.py`:

1. **The arm travels with the value.** `armed_external = bool(external_armed and
   np.isfinite(external_residual) and not generic_departing)`. A channel that reported no usable number
   cannot vouch. Ordinary reactors pass `external_armed=False`/`external_residual=nan`, so this is a no-op
   there.
2. **The fold distinguishes `nan` from `±inf`.** `nan` = no information → keep the generic residual
   (ordinary path, unchanged). `±inf` = the channel is reporting a NON-FINITE residual → POISON the folded
   residual (`r = inf`) so the non-finite guard resets the streak and it can never arm, and NAME the
   invalid value in `worst_label` (diagnostics). Dropping it while the generic channel terminated was the
   defect; poisoning blocks termination through **any** channel, even an armed-and-flat generic one
   (`test_a_nonfinite_external_residual_poisons_even_an_armed_generic_channel`).

**Addendum — the compiled hook could emit `+inf`.** `steady_state_external_residual` (`plasma.pyx`)
guarded the electron MOLES (`ne_now > 0`) but took `log()` of the FRACTION `xe = ne/tot`. A positive
subnormal `ne` divided by an order-one total underflows to `0.0`, so `log(0) = -inf` → the hook returns
`+inf`. Confirmed through the **compiled** hook. Fixed by guarding the fraction after the division
(return `nan`, the channel's "no information" sentinel), so the vocabulary {`nan`, finite} is restored at
the source and a generic `inf` can only ever arrive on the generic side. `rework.md`'s earlier "never
`inf`" claim is corrected above. Post-fix demonstrations: `1e-9` follows the tolerance rule; `nan`,
`+inf`, `-inf` cannot terminate; the value is named; residual/Jacobian recovery is unchanged (Check 7).

### HIGH 2 — a zero denominator makes the ratio criterion undefined (owner's ruling; my earlier fix was forbidden)

The partial round's fix reproduced the historical `edge/0 == +inf` to preserve promotion. The owner
**forbade** that: a zero denominator makes the relative criterion undefined, and it may then grant neither
promotion nor termination; `denominator = 1`, `= epsilon`, or any floor is a category error (a dimensional
rate compared against a dimensionless tolerance — `1.0` was the bug, a smaller number the same bug). The
enlargement ratio is now computed by `ReactionSystem._rate_ratios_or_zero(rates, char_rate)`, which
evaluates `|rates/char_rate|` **only** when `char_rate` is finite and strictly positive and otherwise
returns **zeros** — the criterion abstains, promoting and terminating nothing, laundering no `0/0` NaN.

**The contract already contains the dimensioned absolute criterion** the ruling requires for the
zero-core-flux case: the `_CHAR_RATE_FLOOR` band on `total_char_rate` (the zero-flux promotion block and
the steady-state-inert block), which includes the wall's `get_non_chemical_char_rate()`. So there is **no
policy gap**, and none was invented. Seven-case denominator matrix, all green:

| denominator | ratio criterion | test |
|---|---|---|
| positive (finite) | evaluates `|rates/denom|` (unchanged) | `test_rate_ratio_criterion_evaluates_only_for_a_finite_positive_denominator` |
| zero | abstains → zeros | same (unit) |
| negative | abstains → zeros | same (unit) |
| NaN | abstains → zeros | same (unit) |
| +inf / -inf | abstains → zeros | same (unit) |
| zero-core-flux reactor (non-plasma) | ratio abstains; absolute criterion governs (deterministic `never started`, not a magnitude-dependent promotion) | `test_a_zero_core_flux_reactor_does_not_promote_by_dimensional_comparison` |
| non-zero wall loss + zero core flux (plasma) | ratio abstains; `get_non_chemical_char_rate() > 0` keeps the absolute criterion live | `test_zero_core_flux_with_wall_loss_abstains_on_ratio_but_keeps_an_absolute_criterion` |

### MEDIUM — the report names the criterion that fired

With several steady-state criteria, the readback (`steady_state_residual`) and the success log used
`steady_state_terms[0]` unconditionally, producing statements false on their face (e.g. "below tolerance
1.0000e-30 for 0 consecutive steps" while a different criterion terminated the run). `base.pyx` now
records the term that fired and reports its residual, tolerance and streak; the backstop-stop readback
still falls back to criterion zero. Pinned by `test_steady_state_report_names_the_criterion_that_fired`
(criteria of tolerance `1e-30` and `1e-6`).

### Test findings folded in

- **Source apportionment** (`plasmaWallTest.py`): the pair source is one reaction `Ar -> Ar+ + e-`; the
  test asserted only the electron leg. It now asserts cation production and neutral consumption too
  (each equals the full source at the all-neutral state where the wall loss is zero), and that the
  non-ionisable bath gas stays out of the balance.
- **AST census** (`plasmaWallTest.py`): `has_assert` accepted a vacuous `assert True` (a bare-`Constant`
  test), and the tag regex accepted a bare severity word. Tightened: a constant assertion no longer
  counts, and a severity word must carry its finding number or a colon (`HIGH 1`, `MEDIUM:`), not appear
  alone.
- **Cross-channel census pin** (new, `plasmaWallTest.py`): pins the flux/transport fold sites to source
  and forbids the reintroduction of `char_rate if char_rate > 0.0 else 1.0` or any `ratio_denom` floor,
  and requires the abstaining helper.

### The census detector was fixed, not just extended

The round-109 cross-channel census enumerated the **residual, charge and Jacobian** channels and reported
**6** sites. The detector's *channel taxonomy* omitted the **FLUX / TRANSPORT** channels, so it was blind
to the two folds that combine gas-phase chemistry with a non-chemical rate — one of which carried this
round's HIGH 2. Re-run with the taxonomy completed to include flux/transport, the enumeration is:

1–6. the six round-109 sites (the fold, the arm gate, the combined arm, and the three charge/quasineutrality
   sites), unchanged.
7. **`base.pyx` — the chemical/non-chemical FLUX fold.** `total_char_rate = sqrt(char_rate² +
   non_chemical_char_rate²)` combines the chemistry rate with the transport rate, and the enlargement-ratio
   denominator ranks edge flux against `char_rate`. This is where HIGH 2 lived: the `else 1.0` truthiness
   substitution erased the "no chemistry → undefined ratio" meaning into a dimensional comparison.
   **1 defect before this round; 0 after** (the ratio now abstains; a non-finite flux is refused loudly by
   the existing `base.pyx` non-finite guard, not laundered).
8. **`plasma.pyx` — the chemistry/wall-transport fold.** `_apply_wall_terms(y, V, res)` combines the
   chemistry residual `res` with the wall transport term. A non-finite wall/source term is refused by the
   round-106 source guard; a non-finite `res` propagates as a bad residual the integrator rejects. **0
   defects.**

**Count: 8 cross-channel sites (was 6; +2 flux/transport). Erasure defects: 1 before this round (site 7,
HIGH 2); 0 after. Zero reported as zero.** Why the two were invisible: the detector searched only the
steady-state residual and charge channels; a fold is only "cross-channel" to it if both operands were in
that list, and neither `char_rate`/`non_chemical_char_rate` nor the wall transport term was. The pin test
now holds the flux/transport sites to source so the list cannot silently drift from the code again.

### The bitwise baseline is no longer an untracked file

`verify_zero_wall_bitwise.py` imported `rmgpy.solver.plasma_base_i246`, a 1217-line **untracked** copy of
the pre-change reactor — so the bitwise-equivalence claim was unreproducible from a clean checkout. The
copy was byte-identical to `git show 311818121:rmgpy/solver/plasma.pyx`, i.e. a pure artifact of an
interrupted `build_and_verify_zero_wall.sh` run (that wrapper already git-shows the baseline, builds it,
and cleans up on exit). The stray was removed; the verify script now exits with a clear message pointing to
the wrapper if the transient module is absent, rather than a bare `ImportError`. A frozen private copy is
the wrong instrument precisely because nothing makes it track the ancestor it claims to be.

### Close gate

- Focused suites rebuilt and green: `pytest test/rmgpy/solver/plasmaWallTest.py
  test/rmgpy/solver/steadyStateTest.py test/rmgpy/rmg/inputTest.py -o addopts="" -p no:cacheprovider`
  → **239 passed, 1 skipped** (was 3 failed, 233 passed, 1 skipped on the owner's rebuild of the partial).
- `.so` proven by value: `strings rmgpy/solver/base.*.so | grep _rate_ratios_or_zero` (present); both
  `base` and `plasma` rebuilt.
- **Wall-loss frequency re-measured and unchanged:** at Te = 3000 K, `PlasmaReactor.compute_nu_wall` =
  15.954491 s⁻¹ and the independent closed form = 15.954492 s⁻¹ (agree to 5.6e-8), matching the owner's
  15.954516 s⁻¹ to 1.5e-6 — far inside the ±3% mobility tolerance. The sub-threshold floor is the pure
  quotient `S/nu_wall`; with `nu_wall` unchanged and the source an input, it is unchanged (closed-form
  estimate 2.33e-20, same order as the owner's integrated 2.5543e-20). None of the round-110 edits touch
  `compute_nu_wall`, `_apply_wall_terms`, or the wall coefficient. Evidence:
  `evidence/round110_wall_remeasure.log`.
- **Check 7 (zero-wall bit-for-bit)** still passes: residual, Jacobian and `jacobian_matrix` identical to
  the pre-change build on a wall-less reactor across 7 ionisation degrees, negative control confirms it can
  fail. Evidence: `evidence/round110_zero_wall_bitwise.log`.
- No file under `rmgpy/molecule/`, `rmgpy/kinetics/`, `rmgpy/data/` or the database touched. Nothing
  pushed, merged or rebased.

## Round 111: the infinities were poisoned and NaN was not — two HIGH, one MEDIUM, two LOW

Round 110 closed the non-finite defect for the values it named. Round 111 is the same defect at its
third appearance, on the members round 110 left open — plus a second entry point into the ratio bug
that the round-110 helper install did not route through. `plasma.pyx` is untouched again (the wall
physics did not move — a close condition): only `termination.py` (pure Python) and `base.pyx` changed.

### HIGH 1 — a non-finite residual must be closed as a CLASS, not member by member

The round-110 authorisation matrix ended at `+inf → 0, -inf → 0, NaN → 1 unsafe path`. Arm the generic
channel with `R = 2` over `[1, e]`, then feed two flat generic intervals with `external_residual = NaN`
while the external channel is **armed** — the second interval returned `True`. This is not hypothetical:
`steady_state_external_armed` returns True (electron present, `t·nu_wall ≥ 1`) while
`steady_state_external_residual` returns **NaN** on a subnormal electron-*fraction* underflow — exactly
the case round 110's own addendum **created** when it rerouted that underflow from `+inf` onto the NaN
branch. NaN then took the benign "no information" branch and a flat generic channel fired on an
armed-but-unusable electron channel.

The fix restructures the fold in `update()` so no non-finite value can reach the flat test as a usable
number. The dispatch is `np.isfinite` — whose negation is True for NaN **by construction** — never a
comparison or a maximum (`NaN > r`, `NaN < r`, `NaN == r` are all False, so a comparison-based guard
passes NaN silently). A non-finite external residual poisons the fold with `inf` when the channel is
**in play**: `external_armed` (a live discharge vouching for the electron) OR an actual `±inf` (which
only an active channel can compute — an ordinary reactor emits only the default NaN, so `not np.isnan`
keeps the `±inf` poison unconditional as in round 110). The only surviving benign case is a NaN from a
channel not in play — the ordinary reactor, whose generic-only termination stays byte-for-byte. The
armed-generic test, which covered only the infinities, now carries the whole class (`(inf, unarmed)`,
`(-inf, unarmed)`, `(NaN, armed)`), and a dedicated test reproduces the reviewer's exact sequence.

### HIGH 2 — the surface-species ratio bypassed the helper and still divided by zero

`base.pyx` builds the surface-to-core promotion ratio as `max(|production|, |consumption|) / char_rate`
**directly**, not through the round-110 `_rate_ratios_or_zero`. On a reversible surface reaction with
equal forward and reverse flux the net rates cancel (`char_rate == 0`) while the gross rates are
positive, so the division is `positive / 0`: the production entry **raised `ZeroDivisionError`** (a numpy
divide would instead have promoted the surface species through an undefined criterion). Reproduced end to
end on the production path — a reversible `A ⇌ B` whose species carry identical thermochemistry (K_eq = 1)
so equal moles give equal fluxes, `B` declared a surface species — the pre-fix `base.so` crashes at
`base.pyx:1123`. Routed through `_rate_ratios_or_zero` (finite, positive denominator or abstain — no
`1.0`, epsilon or floor, the owner's ruling as for core/edge), the ratio abstains, nothing is promoted,
and the run completes. The **positive-denominator arm** (a small thermodynamic offset → `char_rate > 0`)
still promotes `B` on the same path, so the fix removed only the undefined case — a guard broken to never
promote is not what shipped.

### The census detector was repaired, not extended (LOW 1)

The round-110 cross-channel census asserted source **substrings** for the folds it already knew about —
so a comment carrying the string satisfied it, a behaviourally-identical rewrite broke it, and, decisively,
it could not see a division site nobody had listed. That is why it missed the surface ratio: the census's
own round finding it could not detect. The repair **enumerates** rather than lists. Every code division by
`char_rate` in `base.pyx` is found structurally via the tokenizer (`base.pyx` is not valid Python, so it
tokenizes rather than `ast.parse`-s; comments and strings are dropped, so prose and dead code cannot
register), and each must fall inside `log_rates` — the only legitimate site, a display line guarded by
`if char_rate == 0.0`. The surface bug lands outside it and fails the census; after the fix the count
outside is zero. This is the round-111 HIGH 2 red-first: the detector fails on the exact unguarded
division at line 1123 before the fix, passes after.

### MEDIUM — `update()` still leaked a numpy.bool_ on the relaxation-time fallback

The round-110 source fix cast the `armed_external` attribute, but the default e-fold fallback computes
`span_ok = (np.log(t_now) - np.log(t_flat_start)) >= 1.0` — a `numpy.float64` comparison yielding a
`numpy.bool_` — and returned `self.armed and span_ok` unchanged, so a genuine termination on that path
returned a `numpy.bool_`, not a Python bool. The return is now cast at the single chokepoint every
return-True path flows through (`return bool(self.armed and span_ok)`). A test fires a termination on the
fallback and asserts `type(v) is bool`.

### The census could not fail itself (LOW 2), and caught a live test

The round-106 test-hygiene census accepted a docstring naming a round/finding tag **OR** an assertion, so
a tagged docstring with a `pass` body satisfied it while asserting nothing. A tag is prose about intent;
only an assertion backs a claim. Every `test_*` in the two solver files must now contain a real
(non-vacuous) assertion, tag or no tag. A tripwire feeds a tagged-`pass` and a bare `assert True` and
confirms both are reported. The tightening immediately caught a **live** assertion-less test —
`test_wall_only_run_does_not_divide_by_zero_chemistry_rate`, which relied on "does not raise" with no
`assert` — now given an explicit assertion (the run completes to its backstop, premise `char_rate == 0`
guarded).

### Close gate

- Focused suites rebuilt and green: `pytest test/rmgpy/solver/plasmaWallTest.py
  test/rmgpy/solver/steadyStateTest.py test/rmgpy/rmg/inputTest.py -o addopts="" -p no:cacheprovider`
  → **244 passed, 1 skipped** (was 239 passed, 1 skipped; +5 tests). Evidence:
  `evidence/round111_reds.log`, `evidence/round111_greens.log`, `evidence/round111_suites.log`.
- Only `base.pyx` (compiled) and `termination.py` (pure Python) changed; `base.so` rebuilt (`make build`,
  exit 0). `.so` proven by value two ways: the surface reproduction flips `ZeroDivisionError` → clean on
  the rebuilt module, and `strings base.*.so | grep _rate_ratios_or_zero` shows the sanctioned helper
  compiled in.
- **Wall-loss frequency re-measured and unchanged:** at Te = 3000 K, `compute_nu_wall` = 15.954491 s⁻¹,
  independent closed form = 15.954492 s⁻¹ (agree to 5.6e-8), matching the owner's 15.954516 s⁻¹ to 1.5e-6;
  the sub-threshold floor is the pure quotient `S/nu_wall`, unchanged (2.33e-20, same order as the owner's
  2.5543e-20). `plasma.pyx` is untouched in the round-111 diff. Evidence: `evidence/round111_remeasure.log`.
- No file under `rmgpy/molecule/`, `rmgpy/kinetics/`, `rmgpy/data/` or the database touched. Nothing
  pushed, merged or rebased.

## Round 112: one gate for every non-finite, at the decision boundary — two HIGH, one LOW

The non-finite defect at its fourth appearance. Rounds 109→110→111 closed it value-by-value; round 112
closes it by CONSTRUCTION at each decision boundary, so a value the enumeration never named is refused
without being enumerated. `plasma.pyx` is byte-identical (verified against both `HEAD` and the round-110
base `ca384f8f3` — the wall physics did not move, a close condition); only `termination.py` (pure Python)
and `base.pyx` (rebuilt) changed.

### HIGH 1 — the steady-state decision: absence made distinct from a supplied non-finite value

Two reproductions of one defect. (a) `termination.py` poisoned a `nan` external residual only when the
channel was `external_armed=True`; with `external_armed=False` a supplied `nan` slipped onto the benign
"no information" path and a flat, armed generic channel fired on it (`False, True` over two flat
intervals). (b) `_slope_analysis` took each live species' log-log slope and **dropped** any that came
out non-finite (`finite = np.isfinite(slope)`, keeping only `finite[j]` entries): with moles
`[1e308, 1e-20]` and a floor below `1e-20` the trace species is live but its mole fraction
`1e-20 / 1e308` underflows to zero, `log(0) = -inf`, `inf - inf = nan`; the finite heavy species then
governed and, once armed and flat, two intervals returned `False, True`.

The root of (a) was that a bare `nan` VALUE could not be told apart from "no external channel" — so the
guard keyed on the arm FLAG. The fix makes **absence a distinct sentinel**: `external_residual` defaults
to `None`, structurally distinct from any supplied non-finite value. When a channel is supplied (any
non-`None` value), a non-finite reading — `nan` or ±inf, armed or not — poisons the folded residual with
`inf` so the single existing guard `if not np.isfinite(r): self.streak = 0; return False` refuses. No
per-case branch on `isnan`/`isinf`/`armed`; the value is never dropped. For (b), a non-finite per-species
slope now returns status `'inf'` (the whole step poisoned, `r = inf`), naming the offending species,
exactly as an appearing or large-negative species already was — closed with `np.isfinite(slope).all()`,
never a per-element comparison a `nan` slips through.

`base.pyx` maps its overloaded plasma hooks onto the new sentinel at the call site: the external channel
is in play only when it has something to say this step — `ss_external_armed or np.isfinite(ss_external_
residual)` — otherwise `None`. This reproduces the pre-round-112 INTEGRATION behaviour byte-for-byte (an
ordinary reactor, and a wall deck with no ionisation source whose electron merely decays, both report an
unarmed `nan` and are absent → the generic channel decides), while the unit-level guarantee is newly
enforced for direct callers. **The first cut of this map keyed presence on a finite `relaxation_time`
and regressed the gamma=0 wall-only deck** (`test_wall_only_deck_integrates_past_t0_on_the_production_
path`): that deck has a finite `relaxation_time` but an unarmed `nan` residual, so the criterion refused
termination, the run over-integrated, and the electron drifted to `-3e-18 mol` and tripped the wall
guard. "Has something to say" is the correct presence signal; the caught regression is recorded because
the discriminator is exactly the kind of adjacent-quantity read this campaign keeps finding.

### HIGH 2 — the enlargement/interrupt/promotion ratios: one finiteness gate they all route through

A finite, positive denominator does not guarantee a finite ratio. `_rate_ratios_or_zero` returned
`np.abs(rates / denominator)` unchecked, so a ±inf rate (an unresolved network-leak rate never covered by
the raw-rate guard at `base.pyx:931`, which checks only `char_rate` and the edge rates) produced an
`inf` ratio, and a finite-but-huge surface gross rate **overflowed** to `inf` even over a finite
denominator. Either non-finite ratio then defeats the gates that read it — `inf > tol_interrupt`
interrupts the run and promotes numerical garbage into the core (`base.pyx:1216`, `:1307`), and `nan`
slips every comparison yet is picked as the `argmax` maximum (surface promotion, `base.pyx:1136`).

Core, edge, network-leak and surface ratios **all route through `_rate_ratios_or_zero`** — the round-111
census proves no site bypasses it — so the single mechanism lives there: over a finite, positive
denominator, if any computed ratio is non-finite the helper stops loudly (`ValueError`), naming the
offending indices and their rates, rather than promote it or launder it to zero. The round-110/111
abstention on a zero/negative/non-finite denominator (return zeros) is preserved, and the existing
non-finite-RATE guard at `:931` still fires first for a `nan` `char_rate` (the helper abstains there).
A new input path into any of the four sites is covered by construction.

### LOW — the assertion census still had four holes

`has_real_assertion` walked the whole function with `ast.walk` and accepted: `assert 1 == 1` (a constant
COMPARE — the round-111 census only rejected a bare `ast.Constant`); an `assert` under `if False:` that
never runs; an `assert` inside an uncalled nested function that `ast.walk` reaches but the body never
executes; and a context manager whose AST merely CONTAINS the substring `raises`/`warns` (e.g. a helper
named `a_thing_that_raises`). The detector now walks reachable statements only — never into a nested
def/class or a statically-dead branch — treats an assert whose test references no runtime value as
vacuous, and recognises `pytest.raises`/`.warns` by the call TARGET, not a substring. The census's own
self-test gained the four reproductions (red before, green after); running the tightened census over the
two solver test files flagged no live test this round.

### Close gate

- Full solver suite rebuilt and green: `pytest test/rmgpy/solver/ test/rmgpy/rmg/inputTest.py -o
  addopts="" -p no:cacheprovider` → **348 passed, 1 skipped** (`test/rmgpy/solver/` alone: **247 passed**,
  0 collection errors, all 9 files collected). The new tests are red on `8dcfcd223` and green after:
  `evidence/round112_reds.log` (3 failed pre-fix), `evidence/round112_greens.log`.
- Only `base.pyx` (compiled) and `termination.py` (pure Python) changed; `base.so` rebuilt (`make build`,
  exit 0). `.so` proven by value: `strings base.*.so | grep 'Non-finite enlargement rate ratio'` shows
  the new gate string compiled in, and the HIGH 2 reproduction flips DID-NOT-RAISE → `ValueError` on the
  rebuilt module.
- **`plasma.pyx` byte-identical** (`git diff HEAD -- rmgpy/solver/plasma.pyx` and `git diff ca384f8f3 --
  …` both clean). **Wall-loss frequency re-measured and unchanged:** `compute_nu_wall` = 15.954491 s⁻¹,
  closed form 15.954492 s⁻¹, sub-threshold floor 2.33e-20 — identical to round 111. Evidence:
  `evidence/round112_remeasure.log`.
- No file under `rmgpy/molecule/`, `rmgpy/kinetics/`, `rmgpy/data/` or the database touched. Nothing
  pushed, merged or rebased.

## Round 113: the adapter and the denominator — the same hole, one layer out — two BLOCKING, one LOW

Round 112 closed both reproductions inside `TerminationSteadyState.update()` (they now give `False, False`)
and routed all four ratio paths through the helper. The spar still did not clear: production goes through
two layers the update()-level tests bypass — the base.pyx **adapter** that maps the reactor's hooks onto
the criterion, and the helper's **zero-denominator abstention**. This round tests THROUGH `base.pyx`, not
by calling `update()`. `plasma.pyx` is byte-identical again.

### BLOCKING 1 — the adapter laundered a source-driven channel's non-finite reading to absence

The round-112 adapter keyed the external channel's presence on `external_armed OR isfinite(residual)`, so
an UNARMED non-finite reading (a `nan` while the electron is unresolvable, or a blown ±inf) was mapped
back to `None` — "no channel" — and a flat generic channel fired on it: the same laundering as round 111,
one layer out. Whether the channel is in play is a **structural** fact about the reactor — does it drive an
electron from an ionisation source — not a property of the value it reported this step. The adapter now
reads that: a reactor with an active source has a channel every step, so ANY non-finite reading it produces
poisons (armed or not); a reactor with no source — an ordinary reactor, or a wall deck whose source-less
electron merely decays — has no such channel and its `nan` is genuine absence (plasma's documented "no such
channel" sentinel), so the generic channel decides alone.

This is the crux the round-112 regression already exposed: plasma's hook returns `nan` for BOTH "no source
channel" and "channel present but unresolvable", and the electron-underflow wall guard that crashes an
over-integrated afterglow lives in the frozen `plasma.pyx`. So the gamma=0 wall-only deck
(`test_wall_only_deck_integrates_past_t0_on_the_production_path`) and the reviewer's poison test feed the
adapter an **identical** value (`nan`, unarmed, finite relaxation time) yet require opposite outcomes —
absence vs poison. The only bit that separates them is source presence, which plasma's overloaded `nan`
cannot carry. The adapter reads `ionisation_source` from the reactor (read-only `getattr`, no `plasma.pyx`
change) to recover it. The truly clean cure — a hook returning a distinct absence sentinel — lives in
`plasma.pyx`, out of scope; the coupling is noted so the owner can rule on it. Tested through the
production `simulate()` path with a `PlasmaReactor` subclass whose hook returns each poison unarmed, driven
to a flat generic channel (red on the round-112 `.so`: the adapter handed the criterion `None`; green
after: it hands the poison, and the run does not reach steady state).

### BLOCKING 2 — `_rate_ratios_or_zero` laundered an infinite numerator to zero

The helper returned zeros whenever the denominator was zero/non-finite, EVEN IF a numerator was infinite:
`core_species_rates=[0]`, `char_rate=0`, `network_leak_rates=[inf]` → `[0]`, and the inert/steady-state or
interrupt path then acted on a zeroed inf. The raw-rate guard did not inspect network rates. A non-finite
**numerator** is a broken integration regardless of the denominator, so the helper now stops loudly first,
before the denominator is consulted; the zero-denominator abstention applies only to FINITE numerators
(the undefined `0/0`, `x/0`). The raw-rate guard was moved ahead of the ratio calls and extended to inspect
`network_leak_rates` too, so a non-finite network rate stops loudly with the detailed message (preserving
the existing core/edge rate-guard messages by construction, since the guard fires before the helper). The
round-112 assertion that blessed the laundering (`ratios([inf,1.0], 0.0) → zeros`) is inverted: it now
raises, while a FINITE numerator over a zero/non-finite denominator still abstains.

### LOW — NumPy error configuration, and four more census holes

The overflow finiteness check is computed under `np.errstate(all='ignore')` and then tested, so the
verdict does not depend on NumPy's global divide/overflow settings. The assertion census now recognises a
`raises`/`warns` context manager by the call TARGET being `pytest.<raises|warns>` (rejecting a bare
`raises(...)`, which need not be pytest's, and an unrelated `fake.raises(...)`), and treats assertions
under `while False:` and in empty loops (`for _ in []:`, `range(0)`) as unreachable, alongside the round-112
`if False:` and uncalled-nested-function cases. Its self-test gained the four reproductions (red before,
green after); the tightened census flags no live test.

### Close gate

- Full suite rebuilt and green: `pytest test/rmgpy/solver/ test/rmgpy/rmg/inputTest.py -o addopts="" -p
  no:cacheprovider` → **349 passed, 1 skipped** (`test/rmgpy/solver/` alone: **248 passed**, 0 collection
  errors). The new/inverted tests are red on the round-112 tip and green after: `evidence/round113_reds.log`
  (2 failed pre-fix), `evidence/round113_greens.log`.
- Only `base.pyx` (compiled) changed this round; `termination.py` untouched. `base.so` rebuilt (`make
  build`, exit 0). `.so` proven by value: `strings base.*.so | grep 'Non-finite enlargement rate:'` and
  `… | grep 'network indices'` both present.
- **`plasma.pyx` byte-identical** (vs `HEAD` and vs `ca384f8f3`). **Wall-loss frequency re-measured and
  unchanged:** `compute_nu_wall` = 15.954491 s⁻¹. Evidence: `evidence/round113_remeasure.log`.
- No file under `rmgpy/molecule/`, `rmgpy/kinetics/`, `rmgpy/data/` or the database touched. Nothing
  pushed, merged or rebased.
