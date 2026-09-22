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

## Files touched

`rmgpy/solver/plasma.pyx`, `rmgpy/solver/base.pyx`, `rmgpy/solver/base.pxd`,
`rmgpy/solver/termination.py`, `rmgpy/rmg/input.py`, `documentation/source/users/rmg/input.rst`,
`test/rmgpy/solver/plasmaWallTest.py`, `test/rmgpy/solver/steadyStateTest.py`, and this
`docs/i246-ambipolar-wall-operator/` directory. No file under `rmgpy/molecule/`, `rmgpy/kinetics/`,
`rmgpy/data/` or the database was touched. Nothing pushed, merged or rebased.
