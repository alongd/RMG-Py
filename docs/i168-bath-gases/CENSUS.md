# I-168 — who depends on RMG's automatic bath gases?

Since commit `5ec0636e` (Connie Gao, 2014-08-04) every RMG job has had four species pushed
into it that no input file asked for:

```python
# Also always add in a few bath gases (since RMG-Java does)
for label, smiles in [("Ar", "[Ar]"), ("He", "[He]"), ("Ne", "[Ne]"), ("N2", "N#N")]:
    molecule = Molecule().from_smiles(smiles)
    spec, is_new = self.reaction_model.make_new_species(molecule, label=label, reactive=False)
    if is_new:
        self.initial_species.append(spec)
```

The commit message, in full, is *"Add set of bath gases (Ar, Ne, He, and N2) to all RMG jobs
just like Java."* Every subsequent touch of these lines (`d74b859f8a`, `14bfbbe163`,
`4a4173993a`, `e61e0d52e5`) is Python-3 or PEP-8 churn. **In eleven years no commit has stated
a technical reason for the line, and no user documentation has ever mentioned that it happens**
(`grep -rni "bath gas\|always add" documentation/source/users/` returns only Arkane's explicit
`bathGas` dicts and unrelated prose).

This document is the census that licensed gating it. It is the deliverable, not a footnote:
the change is four lines and the whole risk was in what else leaned on them.

## Why it matters — what the species actually do

They are created `reactive=False`, so they generate no chemistry. But `initialize()` enlarges
**nonreactive species into the core first**, ahead of the deck's own:

```python
# Add nonreactive species (e.g. bath gases) to core first
# This is necessary so that the PDep algorithm can identify the bath gas
for spec in self.initial_species:
    if not spec.reactive:
        self.reaction_model.enlarge(spec, requires_rms=requires_rms)
```

So they are core species with `index == -1`, which sorts them to the head of every Chemkin
`SPECIES` block, species dictionary, Cantera YAML and saved seed mechanism. Measured on the
5 torr argon deck (`docs/i164-ar5torr/input.py`), before the change:

```
SPECIES
    He        Ne        N2        e-(1)     Ar(2)     Arp(3)     [Li](4)     [Lip](5)
END
```

Three of the eight are species the modeller never mentioned, each at exactly zero mole
fraction. (`Ar` is *not* one of them — the deck declares argon itself, so `is_new` is False and
the injected `Ar` is discarded. `[Li]`/`[Lip]` come from a reaction library and are a separate
ticket.)

## Method

Three independent searches, each reporting hits classified **hard** (would break), **soft**
(behaviour changes silently) or **incidental** (name appears, no dependency):

1. `test/` — every test, fixture and expected-output file.
2. `rmgpy/`, `arkane/`, `scripts/` — every runtime path, `.py`/`.pyx`/`.pxd`.
3. `examples/`, `test/regression/`, `documentation/`, repo root — every RMG input deck,
   asking whether it *references* one of the four labels without *declaring* it.

Plus the git history of the line itself.

False positives were filtered on word boundaries: `Ar` is a substring of `Arrhenius`,
`ArrheniusBM`, `Aromatic`, `Area`; `He` of `Heat`, `Here`; `Ne` of `New`, `Never`, `Net`.

## Result — the injection has real dependants, and none of them are plasma

### Hard (removal breaks them)

| # | Site | What breaks |
|---|---|---|
| 1 | `rmgpy/rmg/pdep.py:865-866` | `bath_gas = [spec for spec in reaction_model.core.species if not spec.reactive]` followed by `assert len(bath_gas) > 0`. A pressure-dependent job whose deck declares no inert of its own is kept alive **solely** by the injection. |
| 2 | `rmgpy/solver/simple.pyx:253-257` and `mbSampled.pyx:217-221` | `calculate_effective_pressure` indexes `y0` (length `num_core_species`) by `species_index[rxn.specific_collider]` (built over core **and** edge). A `(+N2)`-style collider is only in range while N2 is a core species. The equivalent lookup in `residual` (`simple.pyx:411-413`) carries a length guard; this one never got it. |
| 3 | `test/rmgpy/test_data/restartTest/*/filters/species_map.yml` | Both fixtures enumerate Ar/He/Ne/N2 as the first four entries (26 mapped = 22 seed + 4; 9 mapped = 5 seed + 4). `main.py` raises `RuntimeError` for a mapped species missing from the core, so `TestRestartWithFilters` and `TestRestartNoFilters` in `test/rmgpy/rmg/mainTest.py` would both fail. |

### Soft (removal changes results silently)

* **Pressure-dependence physics.** `rmgpy/rmg/pdep.py:868-871` weights *every* nonreactive core
  species equally — `self.bath_gas[spec] = 1.0 / len(bath_gas)`, carrying the source comment
  `# is this really the only/best way to weight them?`. Consumed by
  `rmgpy/pdep/configuration.pyx:190-201` to build σ, ε and molecular weight. So **today every
  pdep network in RMG is computed against an equimolar Ar/He/Ne/N2 mixture**, and a deck that
  declares its own N2 bath gas gets that N2 diluted to 25%. Removing the injection would move
  the collision frequencies, and hence the rate coefficients, of every existing pressure-
  dependent model. Arguably a correctness fix; unarguably a silent numerical change.
* **Output artifacts.** Element lists (`rmgpy/rmg/model.py:168-186` → `chemkin.pyx:2523-2535`,
  `yaml_cantera2.py:331-350`) lose He/Ne/Ar; third-body efficiency tokens for the four vanish
  from `chem.inp`, Cantera YAML and RMS YAML; saved seed mechanisms shrink.
* **CI regression baselines.** All 11 comparisons in `.github/workflows/CI.yml:284-385` diff
  against a baseline generated on `main`; `scripts/checkModels.py:126,133-143` errors on any
  species present in one model and not the other. `aromatics`, `fragment` and
  `liquid_oxidation` declare no inert at all, so their cores would shrink by four.

### Incidental

~13 runtime sites and ~30 test modules name the four labels without depending on the
injection — `chemkin.pyx:993,1446` re-derives inertness by isomorphism when *reading* a file;
`rmgpy/rmg/model.py:117`, `tools/mergemodels.py:172`, `tools/observablesregression.py:287` are
offline label filters; `data/auto_database.py:126-128` skips nonreactive species *and* runs
~100 lines before the injection; `tools/loader.py`, `scripts/thermoEstimator.py` and
`tools/isotopes.py` never call `initialize()` at all.

### Input decks: a clean negative result

**Zero RMG input decks in the repo rely on the injection.** All 36 decks that name one of the
four labels declare it themselves with `species(label=..., reactive=False, structure=...)`. The
only references-without-declaration are three prose snippets in `documentation/source/users/rmg/
input.rst` (lines 415-428, 475-497) and `surfaces.rst` (lines 67-81), all naming `N2` in a
`simpleReactor`/`surfaceReactor` example — none of them plasma, so all three remain correct.

### Empty classes

* **Plasma.** `rmgpy/solver/plasma.pyx` (925 lines) contains **zero** matches for
  `Ar|He|Ne|N2`, `reactive`, `bath`, `inert`, `collider` or `efficienc`. `PlasmaReactor` does no
  third-body or collider handling whatsoever, and its `initial_mole_fractions` validation
  requires every key to be a core species — the injected gases are never keys.
* **Transport.** `rmgpy/transport.py` has no bath-gas or collider lookup at all.
* **Label-keyed efficiency lookup.** Every efficiency resolution in the tree is
  isomorphism-based (`rmgpy/kinetics/model.pyx:405-423,518-530`). There is no
  `efficiencies["Ar"]` that could `KeyError`.

## Conclusion, and the ordering fact that makes the gate safe

The injection is load-bearing for ordinary gas-phase RMG and inert for plasma RMG. So it is
**gated on the reactor type, not removed** — `RMG.add_default_bath_gases()` returns early when
`RMG.uses_plasma_reactor()` is true.

The one way that gate could have opened a hole is hard dependant #1: a plasma deck that also
enables `pressureDependence` would reach `pdep.py:866` with no inert in its core. It cannot,
because `check_input()` runs at `main.py:577` — **before** the injection at `main.py:723` — and
already asserts:

```python
assert any([not s.reactive for s in reaction_system.initial_mole_fractions.keys()]), \
    "Pressure Dependence calculations require at least one inert (nonreacting) species for the bath gas."
```

Those keys are resolved from `rmgpy/rmg/input.py`'s own `species_dict`, populated only by
`species(...)` directives. So the guard has never been satisfiable by the injected gases: any
pdep deck must already declare its own inert, which then reaches the core and satisfies
`pdep.py:866` on its own. The two guards agree, and the gate sits between them without touching
either.

Residual, narrow: a plasma **and** pdep **and** seed-restart job could in principle reach
`pdep.py:866` on a path `check_input` did not cover. No such deck exists in the repo.

## Findings recorded for other tickets, not fixed here

1. `rmgpy/rmg/main.py:377-379` — the pdep inert assert sits **outside** the
   `for index, reaction_system in enumerate(...)` loop that precedes it, so it only ever checks
   the *last* reaction system. A multi-reactor deck whose first reactor lacks an inert passes.
2. `rmgpy/rmg/pdep.py:868-871` — the equal 1/N weighting of the bath gas is, as its own comment
   asks, almost certainly not the best way to weight them. Combined with the injection it means
   no RMG pdep calculation has ever used the bath gas the modeller specified.
3. `rmgpy/solver/simple.pyx:253-257` — missing the length guard its sibling at
   `simple.pyx:411-413` has. Reasoned from `neq`/`species_index` sizing, not observed firing.
4. `rmgpy/tools/observablesregression.py:287` — `inert_list` contains `'[N#N]'`, which is not
   valid SMILES.
5. `plasmaReactor` has no user documentation at all (`grep -rn plasmaReactor documentation/`
   returns nothing), so there is no place to document this behaviour for plasma users yet.
