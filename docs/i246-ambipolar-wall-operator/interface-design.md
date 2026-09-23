# The M8-A wall-flux / wall-energy interface

The successor power-closure workstream (M8-B) is required by contract to consume "M8-A's stable
wall-flux and wall-energy interface". This is that interface: what an accepted-state consumer reads,
how it is latched, what each quantity promises, and what it deliberately does not.

## What a consumer reads

All fields are attributes on `PlasmaReactor`, written only by `_latch_wall_diagnostics` and readable
after `initialize_model` and after every accepted `advance`:

| field | quantity | unit | honesty |
|-------|----------|------|---------|
| `wall_flux` | per-core-species wall particle rate, signed (charged loss < 0, recycled return > 0) | mol/s | **robust** |
| `nu_wall_latched` | common ambipolar loss frequency at the accepted state | s⁻¹ | **robust** |
| `wall_diagnostics_time` | the accepted time the latch is for (staleness guard) | s | — |
| `wall_electron_energy_flux` | `2·k_B·T_e · (electron loss rate)` | W | **robust** (from T_e) |
| `wall_neutralization_energy_flux` | `Σ (H_ion − H_neutral)·(ion loss rate)` | W | conditional (needs ion+neutral thermo) |
| `wall_ion_energy_flux` | ion directed/sheath energy, `e·V_sheath·(ion loss rate)` | W | **declared-absent** |
| `wall_energy_availability` | `field → 'available' \| 'unavailable' \| 'declared-absent'` | — | the map that makes the above machine-readable |

The availability map is the addition the reviewer asked for: a consumer summing the energy fluxes
must distinguish a term that is **declared-absent by design** (the ion directed energy) from one that
merely **could not be computed** (the neutralisation energy when thermo is missing) from one that is
**available** — in code, without reading this document. A bare NaN cannot carry that distinction; the
map can. The power-closure verifier (which must show "power deposited = volume + wall + flow +
stored-energy changes within tolerance") reads `wall_energy_availability['wall_ion_energy_flux'] ==
'declared-absent'` and accounts for that term explicitly rather than discovering a NaN and guessing.

## How it is latched (no Newton-trial leak, by construction)

`wall_loss_rates` remains **residual scratch**: `_apply_wall_terms` overwrites it on every residual
evaluation, including rejected Newton trials. It is *not* the interface, and a consumer must not read
it.

The interface fields are written *only* by `_latch_wall_diagnostics(y, V, t)`, which recomputes the
wall terms cleanly at an accepted `y` and is called from exactly three places, each an accepted-state
gate outside the residual:

- `initialize_model`, on `self.y0` (the initial composition is an accepted state), at t = 0;
- `step`, on `self.y`, **after** `check_wall_support` has passed;
- `advance`, on `self.y`, **after** `check_wall_support` has passed.

`step` is the one that matters in production: `ReactionSystem.simulate` drives the reactor through
`step`, never `advance`. Round 79 wired the latch into `advance` only, so a real simulation published
diagnostics frozen at the t=0 initialisation latch — a guarantee that never ran, worse than none
because it looked like one (round-83 HIGH 1). Both entries now latch, at the same accepted-state gate.

Because every call site is outside the residual — the same accepted-state gate `check_wall_support`
already sits behind — no rejected trial state can reach the interface. This is a structural guarantee,
not a discipline: there is no code path from a Newton trial to `wall_flux`.
(`test_wall_diagnostics_latched_only_at_accepted_states` pins it: a wild residual moves the scratch
and leaves the latch untouched.)

## The wall-energy quantity, and why the ion term is a NaN sentinel not a number

The total energy carried to the wall per electron–ion pair lost is `E_T = E_e + E_neutralisation +
E_ion`:

- **E_e** (electron thermal): the flux-averaged energy of electrons escaping over the sheath is
  `2·k_B·T_e` per electron. The reactor holds T_e, so this is always available and reported.
- **E_neutralisation** (chemical): the formation-enthalpy drop `H_ion − H_neutral` released when the
  ion is neutralised at the surface, from thermo the recycle map already resolves. It is owed for
  **every ion lost to the wall**, not only the recycled fraction: `wall_recycling` (γ) is a
  mass-return fraction governing where the neutral goes, and a fully-pumped ion (γ=0) still recombines
  with an electron at the surface and deposits its enthalpy there. The flux is therefore
  `Σ (H_ion − H_neutral)·(ion loss rate)`, independent of γ — an earlier version multiplied it by γ,
  so a pumping wall reported 0 W and, worse, labelled it `available` (round-88 HIGH 4). Reported when
  both thermo are present and finite; otherwise `unavailable`.
- **E_ion** (directed/sheath): the ion falls through the sheath and arrives with `e·V_sheath`. This
  needs a sheath potential — and a sheath model is an explicit contract non-goal (§ Non-goals; § 10
  of the binding contract names "sheath or plasma-potential assumptions" as M8-B's charter). It is
  reported as a NaN with availability `declared-absent`.

Fabricating `V_sheath ≈ (T_e/2)·ln(M_i/2πm_e)` here was rejected for three reasons: it is out of
scope; § 10 forbids "fit[ting] absorbed power, wall loss and Te simultaneously to one observed
electron density" as non-identifiable, and a fabricated `V_sheath` becomes exactly such a parameter,
leaving M8-B unable to tell this operator's estimate from its own closure; and the campaign's
recurring defect is the plausible wrong answer no check catches — a placeholder V_sheath propagates
silently and looks like a measurement, where a NaN propagates loudly and looks like what it is.

## What the interface cannot promise

- **n_e / the ionisation fraction.** At the self-sustaining root the electron density is marginally
  stable (particle balance fixes T_e, not the density); the interface reports *rates given the
  accepted state*, and the density itself is pinned only by the power balance M8-B closes. See the
  section-14 answer in `rework.md`.
- **The ion directed energy / sheath potential** — declared-absent, above.
- **Validity outside the low-α, electropositive, unmagnetised, diffusion-limited regime** — already
  documented in `input.rst`.
