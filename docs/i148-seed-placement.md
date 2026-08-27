# I-148 — the electron-placement declaration must survive the seed round trip

The fourth audit's single blocking failure (`docs/i123d-audit.md` §9), promoted to its own
ticket and repaired here. Evidence labels as in the audits: **[R]** code, **[D]** database,
**[M]** measured (command + output in `docs/i148-seed-placement/evidence/` or inline),
**[I]** inference with basis.

## 1. The defect, restated

A reaction declares its electron placement **per side** — `(reactant_count, product_count)` —
because a net scalar cannot distinguish an electron consumed from one produced, and the
incident-order information lives on the reactant side. The declaration is stored in
`FAMILY_ELECTRON_PLACEMENT` **[R]** `rmgpy/electron_placement.py`, keyed on the reaction's
owner label.

Every RMG run writes a `seed/` directory and a runnable `restart_from_seed.py` beside its
output, unconditionally **[R]** `rmgpy/rmg/main.py` (`make_seed_mech`,
`initialize_seed_mech`) — a third serialisation of the model, alongside Chemkin and Cantera.
Writing the seed renames the owning library to the fixed string `seed` ("as these names
should not be dynamic"); on restart the directory is copied and loaded under the fixed label
`restart` **[R]** `rmgpy/rmg/main.py` (restart block in `load_input`ed execution path).
Neither fixed label is a registry key, so before this ticket every plasma reaction that
round-tripped through a seed fell back to the net-derived placement rule
**[R]** `rmgpy/electron_balance.py` (`get_electron_placement_counts`, final lines).

### 1.1 One correction to the audit's account

§9.3 says the original label "survives only as free text in `longDesc`, which nothing
parses". That is wrong at the mechanism level, and the correction is what makes the repair
small: `KineticsLibrary.get_library_reactions` **[R]**
`rmgpy/data/kinetics/library.py` DOES parse the seed writer's
`Originally from reaction library: <label>` line back into `LibraryReaction.library` for
auto-generated libraries — and `add_seed_mechanism_to_core` **[R]** `rmgpy/rmg/model.py`
already consumes that parsed field to load the original library on restart. What nothing did
was consult it for PLACEMENT: all three placement-owner lookups read only
`Reaction.family`, which that same reader overwrites with the container label. Measured
**[M]** on the reloaded seed objects, pre-fix:

    Li => [Lip]    type=LibraryReaction family='seed' library='PlasmaElectronImpactIonization'
    [Lip] => Li    type=LibraryReaction family='seed' library='PlasmaRadiativeRecombination'

So the owner already travels with the reaction, in a parsed field; the defect was that the
placement lookups asked the container instead of the reaction.

## 2. The asymmetry, measured before the fix

**[M]** `python docs/i123d-audit/probe_seed_roundtrip.py <deck>` on this tree
(`evidence/probe-roundtrip-BEFORE.log`), deck = `docs/i123-integration/input.py`, exit 0:

| channel | canonical | via seed (before) | via seed (after) |
|---|---|---|---|
| electron-impact ionisation | (1, 2) | **(0, 1)** — declaration lost | **(1, 2)** — preserved |
| radiative recombination | (1, 0) | (1, 0) — coincidence of the net rule | (1, 0) — by declaration |

Identity followed the placement: before, `is_isomorphic(either_direction=True)` was `False`
for the ionisation pair (seed copy ≠ canonical ⇒ added twice) and `True` for the
recombination pair; after, both pairs compare equal.

### 2.1 The two consequences, reproduced before the fix **[M]**

`python rmg.py <deck>/restart_from_seed.py` (the file RMG itself generated, unedited):
core reached **7 species / 3 reactions** — the ionisation channel twice, marked duplicate —
then died at its first Chemkin save
(`MechanismWriterError: ... "Li(2)=>[Lip](3)+e-(1)" has no electron among its reactants`),
exit 1 (`evidence/restart-before-withlibs.log`).

With `reactionLibraries=[]` (seed the only chemistry source): core **7 species /
2 reactions** — the right model — and the save still failed identically, exit 1
(`evidence/restart-before-seedonly.log`). The seed alone is sufficient; the doubling and the
refused save are two consequences of one root cause.

## 3. The repair: the owner travels with the reaction

**The forbidden fix, refused**: adding `seed`/`restart` as registry keys would hand ONE
placement to every reaction that ever passed through a seed, regardless of what it was, and
would silence exactly the failure that made this defect findable. The registry stays a
closed list of real owners.

**What was done instead** — the reaction's own attributions are the owner source:

- `get_placement_owner(reaction)` **[R]** `rmgpy/electron_balance.py` — resolves the
  declared owner from `Reaction.family` first, then from the `Reaction.library` provenance
  the seed round trip preserves. Family-first precedence is monotonic over the old
  behaviour: every reaction that resolved before resolves identically, and the only
  reactions that newly resolve are those whose `family` names no declaration at all.
- `get_placement_declaration` (export + identity) and `resolve_electron_placement` (the
  plasma reactor) both read the owner through it, so the reactor and the writers cannot
  disagree about whom a reaction belongs to.
- The declaration itself still lives in `FAMILY_ELECTRON_PLACEMENT`, once. What travels
  with the reaction is its OWNER IDENTITY — which the seed writer already emits per entry
  and the seed reader already parses back. No serialisation format change, no new field on
  `Reaction`, no second source of truth for what a placement is.

Why this over a per-reaction placement field: a field the seed writer emits and the reader
restores would (a) change the library entry format both repos share, (b) have to be copied
through every path that clones a reaction, and (c) create a second statement of the
placement that can drift from the registry's — the exact two-converters-drift failure I-126
was. The owner-identity route uses a provenance channel that already exists, is already
parsed, and is already trusted by the restart machinery itself (`add_seed_mechanism_to_core`
loads the original library from `rxn.library`).

**The loud tripwire**: `make_seed_mech` now warns per entry, at seed-WRITE time, when a
reaction whose owner carries a declaration is about to be written without provenance the
reader parses (`warn_if_seed_loses_placement` **[R]** `rmgpy/rmg/main.py`;
`seed_placement_survives` **[R]** `rmgpy/data/kinetics/library.py` mirrors the reader's
branches and must be kept in step with them). A representation the provenance lines do not
cover announces itself in the run that writes the seed, not the run that fails to restart
from it. Demonstrated firing on the real writer **[M]**
`python docs/i148-seed-placement/probe_seed_write_warning.py` → P1 (provenance present)
silent, P2 (provenance stripped) fires `SEED LOSES ELECTRON PLACEMENT`, exit 0.

## 4. Proof on the artifact **[M]**

Deck original run on this tree: core 7 species / 2 reactions; seed written; exit 1 at the
pre-existing, out-of-scope Chemkin→Cantera translation step (`cantera_from_ck`, ck2yaml has
no `TDEP/` support — audit §10.3/§10.5, ecosystem converter, separate ticket).

After the fix, `python rmg.py <deck>/restart_from_seed.py`, unedited:

| | core | save | `chem.inp` vs original | exit |
|---|---|---|---|---|
| restart, libraries declared | 7 species / 2 reactions, both channels "not new" | **saved** | **byte-identical** | 1 — same `cantera_from_ck` step as the original run, nothing else |
| restart, `reactionLibraries=[]` | 7 species / 2 reactions | **saved** | **byte-identical** | 1 — same |

The ionisation equation in the restarted model reads
`Li(2)+e-(1)=>[Lip](3)+e-(1)+e-(1)`, second order, as the Voronov coefficient requires.
Exit 0 is unreachable for this deck by the §10.3 defect alone, which predates and is
untouched by this ticket; the restart now terminates in the identical state as the original
run, with an identical mechanism.

## 5. The negative control **[M]**

`examples/rmg/minimal`, same tree, same database: original run 26 species / 66 reactions,
exit 0; restart from its own generated script 26 species / 66 reactions, exit 0, zero
placement warnings. The restart output is NOT byte-identical to the original run's — the
restart machinery renumbers species compactly (e.g. `[CH3](4)` → `[CH3](2)`) and one rate
coefficient reprints as `7.772030e+08` vs `7.772027e+08` — but that is a property of restart
itself, present on the base: the SAME restart run under the base code (fix stashed) produces
`chem.inp`, `species_dictionary.txt` and `chem_edge.inp` **byte-identical to the fixed
tree's**. This change moves a neutral mechanism by exactly zero bytes.

## 6. The edge seed

The shipped deck's edge seed is empty (its whole edge is promoted to the core). A variant
deck (`toleranceMoveToCore=1e6`, otherwise identical) could not hold the charged channels in
the edge either: the plasma solver promotes the cation at t=0 — *"species [Lip](3) was added
to model core to avoid singularity"* **[M]** (`evidence/edge-deck-promotion.log`) — so **no
shipped plasma deck can produce a charged edge seed**. The charged edge-seed round trip was
therefore driven at the object level through the REAL writer and reader **[M]**
`python docs/i148-seed-placement/probe_edge_seed_roundtrip.py` → exit 0
(`evidence/probe-edge-seed-roundtrip.log`): both channels placed in the edge of a hand-built
model, `RMG.make_seed_mech` writes `seed_edge` (its own loop, its own provenance
`try/except`), `KineticsDatabase.load_libraries` + `get_library_reactions` read it back;
family is the container label `seed_edge`, library is the original owner, placements are
(1, 2) and (1, 0), rate-law CLASSES asserted (`VoronovEIArrhenius`, `BadnellRRArrhenius`),
and the write-time warning stayed silent. The neutral minimal example's large non-empty edge
seed additionally round-tripped through the whole-run path in §5.

Not reached end-to-end: a whole RMG run whose final edge holds a charged reaction (no
shipped configuration can produce one), and the filter tensors beyond what the two restarts
consumed. One structural fact found on the way **[R]** `rmgpy/rmg/main.py` (restart block):
on restart the edge seed library is loaded into the kinetics database but the line that
would append it to the model's reaction libraries is commented out
(`# self.reaction_libraries.append(("restart_edge", False))`), so edge-seed reactions are
carried by the restart only as a loaded library, not re-added to the model — the edge seed
is a shallower artifact than it looks.

## 7. Tests

`test/rmgpy/i148SeedPlacementTest.py` — 17 tests: owner recovery (seed and restart shapes,
family precedence, undeclared owners bit-identical to the net rule), identity (the audit's
S4 pair plus the I-134 non-collapse guard `(1,2)ᵀ = (2,1) ≠ (1,0)`), export placement, and
the seed-write warning (fires / silent-with-provenance / silent-for-neutral). The
behavioural red for the recovery tests, pre-fix, is the audit probe itself
(`probe-roundtrip-BEFORE.log`); the file also failed red at collection before the fix
(missing symbols).

Suite, serially, same command (`python -m pytest -m "not functional and not database"`)
**[M]**:

| | failed | passed | skipped | deselected |
|---|---|---|---|---|
| base (ce93e6c12, code = audit base f17bec729) | 2 | 3028 | 49 | 143 |
| with this ticket's commits | 2 | **3045** | 49 | 143 |

The delta is exactly the 17 new tests; the two pre-existing reds on the base
(`electronPlacementTest::test_declaration_registry_is_explicit_and_closed`, whose closed
three-entry pin predates I-119's fourth registry entry, and
`preflightDeckFamilyExclusionTest[docs/i102-quarantine/input.py]`) are unchanged and belong
to their own tickets, not this one.
