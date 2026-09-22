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

## Files touched

`rmgpy/solver/plasma.pyx`, `rmgpy/solver/base.pyx`, `rmgpy/solver/base.pxd`, `rmgpy/rmg/input.py`,
`documentation/source/users/rmg/input.rst`, `test/rmgpy/solver/plasmaWallTest.py`, and this
`docs/i246-ambipolar-wall-operator/` directory. No file under `rmgpy/molecule/`, `rmgpy/kinetics/`,
`rmgpy/data/` or the database was touched. Nothing pushed, merged or rebased.
