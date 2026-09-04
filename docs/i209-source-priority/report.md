# i209 — "the sourced value wins": a guarantee, or an accident of ordering?

Read-only investigation on branch `i209-priority` (`plasma` HEAD `4910d8c52`). No `rmgpy/`
behaviour was changed. Deliverables: this report, a reproduction script
(`docs/i209-source-priority/repro.py`), and a characterization test
(`test/rmgpy/i209SourcePriorityTest.py`, 4 tests, passes on unmodified `plasma`).

## Verdict (Q3), verbatim

**The estimate can win, under these conditions:** a reaction that a kinetics library covers is
*also* carried — with an estimated (or otherwise different) rate — by a **seed mechanism**, and
the seed is loaded in the same run as the library. RMG registers every seed mechanism into the
core (main.py:790-791) **before** it registers any reaction library into the edge
(main.py:795-796), and the duplicate check keeps whichever reaction was registered first while
ignoring kinetics and provenance. The seed's estimate is therefore the incumbent and the
library's sourced rate is the discarded newcomer — the sourced value never becomes operative.
This is not a corner case a user has to construct by hand: RMG's own `--restart` feature loads
the previous run's saved core **as a seed mechanism** (main.py:734,
`self.seed_mechanisms.append("restart")`), so any restart of a run in which reaction R was a
family estimate, after a library covering R has been added, silently keeps the estimate. The
ordering is hard-coded (seeds always precede libraries), so no input configuration avoids it;
only the *absence* of a seed carrying an overlapping estimate avoids it.

## The mechanism (Q1) — confirmed, from code and from a run

`are_identical_species_references` (rmgpy/rmg/model.py:2356) decides "same reaction" from **only**
two things: the reactant/product species references (in either direction, plus `specific_collider`)
and a per-side electron-placement count from `get_electron_placement_counts`
(rmgpy/electron_balance.py). For every neutral reaction the placement count is `(0, 0)` on both
sides, so the verdict reduces to a pure heavy-species reference comparison — kinetics, provenance,
and measured-vs-estimated status are never read. The function's own docstring states the
consequence: *"the model then silently keeps whichever was offered first."*

The keep-the-incumbent behaviour lives in `make_new_reaction` (model.py:536): after resolving
species it calls `check_for_existing_reaction`, and

```python
        if check_existing:
            found, rxn = self.check_for_existing_reaction(forward)
            if found:
                return rxn, False          # model.py:599-602 — returns the INCUMBENT, drops `forward`
```

**Which one the model keeps, and why.** In a normal run the reaction library is registered
(main.py:795-796) before any family reaction is generated (the enlarge loop, main.py:867-872), so
the *library* reaction is the incumbent and the *family estimate* is the newcomer that gets
dropped — the model keeps the **sourced library value**. This is decided purely by registration
order, not by any preference for a sourced rate: `check_for_existing_reaction` returns the first
reaction registered over those species regardless of its kinetics or owner.

### Reproduction (exact code and output)

`docs/i209-source-priority/repro.py`, run in `rmg_env` against the committed
`test/rmgpy/test_data/testing_database` (families `H_Abstraction`; libraries `GRI-Mech3.0`,
`ethane-oxidation`). It builds two reactions over the same heavy species
`[H] + C=C[CH2]C <=> C=C=CC + [H][H]`: a library reaction with a "sourced" Arrhenius
(`A = 1.0e12 cm^3/mol/s` -> `1.0e6` SI) and a reaction with an "estimate" Arrhenius
(`A = 9.9e9 cm^3/mol/s` -> `9900` SI). CASE1 registers the library first; CASE2 registers a
seed-borne estimate first.

```
families: ['H_Abstraction']
libraries: ['GRI-Mech3.0', 'ethane-oxidation']

CASE1 are_identical: True
CASE1 found: True kept is library: True kept kinetics A: 1000000.0000000001

CASE2 are_identical: True
CASE2 found: True kept is seed(estimate): True kept kinetics A: 9900.000000000002 SOURCED A: 1000000.0000000001 ESTIMATE A: 9900.000000000002
```

CASE1: the sourced library value (`A = 1.0e6` SI) survives. CASE2: with the estimate registered
first, `check_for_existing_reaction` returns the **estimate** (`A = 9900` SI) and the sourced
library value is discarded — the same predicate, opposite outcome, decided only by order.

The characterization test `test/rmgpy/i209SourcePriorityTest.py` pins this, including an
end-to-end path through the production entry point `make_new_reaction` where the second reaction
is built from **fresh** species objects (so the collapse is resolved through isomorphism
de-duplication exactly as in a real run, not by sharing Python references):
`test_make_new_reaction_returns_the_incumbent_estimate` asserts `is_new is False` and that the
returned reaction carries the estimate's kinetics.

## Ordering (Q2) — file:line that settles it

`rmgpy/rmg/main.py`, `RMG.initialize`:

- **main.py:790-791** — `for seed_mechanism in self.seed_mechanisms: ... add_seed_mechanism_to_core(...)`
  registers seed mechanisms into the **core**.
- **main.py:795-796** — `for library, option in self.reaction_libraries: ... add_reaction_library_to_edge(...)`
  registers reaction libraries into the **edge**.
- **main.py:867-872** — the `enlarge(spec)` calls that run family generation happen **after** both
  of the above.

So the registration order is fixed: **seed mechanisms -> reaction libraries -> family reactions.**
Each `add_*` method calls `make_new_reaction`, which registers the reaction (model.py:621) if new
and returns the incumbent if not. This settles both halves:

- Library **before** family (795-796 before 867-872) -> a family estimate always loses to a library
  that covers the same reaction. This is the prior review's claim, and it holds.
- Seed **before** library (790-791 before 795-796) -> a seed's rate always wins over a library that
  covers the same reaction. This is the inversion.

**Is there any supported configuration where the ordering differs?** No. The two loops are
unconditional and hard-coded in `initialize`; an input file populates `seed_mechanisms` and
`reaction_libraries` as separate lists but cannot interleave or reorder them, and seeds are always
processed first. The "sourced wins" property for library-vs-family is therefore genuinely
guaranteed by ordering — but so is the "seed estimate wins over library" inversion; it is not a
configuration a user opts into, it is baked into the fixed sequence and triggers whenever a seed
happens to carry an overlapping estimate.

## The two attack candidates (Q3), each probed

**(a) Seed carrying a family-estimated reaction, loaded into core ahead of a library that also
covers it — DOES invert.** Confirmed by CASE2 above and by the end-to-end `make_new_reaction`
characterization test. The realistic trigger is iterative refinement or `--restart`: run 1 has no
library for reaction R, so R gets a family estimate; the model is saved; run 2 adds a library with
the measured R **and** reuses run 1's model as a seed (or is a `--restart`). The seed's estimate is
the incumbent; the library's measured R is discarded. Seed reactions are `LibraryReaction`
objects, so their *type* is "library", but the *numbers* they carry are whatever the prior run
wrote — an estimate. The bug is a value inversion, not a type confusion.

**(b) A family loaded without the kinetics library covering part of its chemistry — does NOT
invert the guarantee.** Here the family estimate does become the operative rate, but nothing is
displaced: there is no sourced value present in the run to lose to it. This is the normal,
intended "no library => use the estimate" behaviour, and it is the user's explicit choice not to
load the library. It is reported here as *tried and does not demonstrate an inversion* — the
estimate is operative by absence, not by silently beating a sourced value that was also loaded.

## Loudness (Q4) — asymmetric, and never names the substitution

What a user sees depends on which reaction is the discarded newcomer:

- **Family reaction discarded (a library or seed was the incumbent) — SILENT.**
  `process_new_reactions` (model.py:972-975): `if not is_new: ... continue`. No log line at any
  level. This is the normal-ordering case, so the ordinary "library beats estimate" collapse is
  completely invisible.
- **Library or seed reaction discarded (something was already registered) — INFO only.**
  `add_seed_mechanism_to_core` (model.py:1760-1761) and `add_reaction_library_to_edge`
  (model.py:1888-1889) both do `if not isNew: logging.info("This library reaction was not new: {0}".format(rxn))`.
  This fires in the inversion case (the sourced library reaction is the discarded newcomer), but it
  is an INFO-level, non-alarming message that says only "not new" — it does **not** say that a
  measured rate was dropped in favour of an estimate, does not name the surviving reaction, and
  does not compare kinetics. In a run whose log is dominated by INFO chatter it is effectively
  invisible.

In neither case does the Chemkin output or the final model report flag the substitution: the
surviving reaction is written with its own (incumbent) kinetics and provenance comment, and the
discarded reaction and its rate simply never appear. There is no warning and no diff. A user would
have to already suspect the collision and grep the log for "not new" to find any trace, and even
then the message would not tell them a sourced value had been overridden.

## Blast radius (Q5)

**Method.** Counting was done over the plasma-lineage database `/home/alon/Code/RMG-database-plasma/input/kinetics`
(the checkout this worktree's `rmgrc` points at), from the source files only — no generation run.
Library reaction population: one reaction = one `^entry(` at column 0 across every `reactions.py`
under `libraries/` (`find libraries -name reactions.py -exec grep -h '^entry(' {} \; | wc -l`;
verified against `grep -c` per file). Families: `ls -d families/*/`. Plasma fraction: sum of
`^entry(` counts over the libraries whose names denote plasma/ion/electron/alkali/Li chemistry.

**Counts.** 190 `reactions.py` files across the library tree (71 standalone + 119 nested inside 10
container libraries), holding **35,216 library reaction entries**, against **138 reaction
families**. The largest single libraries are CurranPentane (3069), CF2BrCl (2455), CH3Cl (2155),
JetSurF2.0 (2087), Chung_solvation_corrections (1784), JetSurF1.0 (1438), NOx2018 (1321). The
genuinely plasma/ion/Li libraries total **207 entries — ~0.6%** of the population; the plasma/Li
families are ~16 of 138 (~12%). So the file-measurable ceiling on the blast radius is **35,216
library reactions x 138 families**, of which ~99.4% is mainstream hydrocarbon / combustion / NOx /
halogen / sulfur / surface chemistry.

**Interpretation.** The collision shape is **general, not Lithium/plasma-specific.** Mainstream
combustion libraries overlap directly with mainstream families (`H_Abstraction`,
`Disproportionation`, `R_Recombination`, `R_Addition_MultipleBond`, `intra_H_migration`) over
exactly the small-radical chemistry those libraries are packed with — the collision was merely
*noticed* in plasma work. There is a second face of the same bug, **library-vs-library**: 50 of the
190 library files contain an H+O2-type branching reaction, 50 contain a CO+OH reaction, 90
reference CH3; canonical small-molecule reactions are the shared backbone of essentially every
combustion mechanism, so two loaded libraries (e.g. GRI-Mech, FFCM1, JetSurF, CurranPentane) will
routinely carry the same reaction with different rates, and load order alone decides which
survives.

**What this count could NOT reach.** The *actual* number of library<->family collisions is a
run-time quantity: family scope is decided by recipe application + group-tree matching on real
molecular graphs during generation, and species-reference identity (and the `(0,0)` electron-count
equality) is established only while a model is built. The 35,216 x 138 figure is a ceiling on the
population that can collide, not a count of pairs that do; the true overlap cannot be measured from
the source files. Library-vs-library identity likewise cannot be confirmed from files (each library
defines its own species objects) — the 50/50/90 figures are label/formula proxies, not confirmed
identity matches.

## What the probe could NOT reach

- **No full RMG generation run was performed.** The evidence is (i) code reading of the identity
  predicate, the duplicate path, and the initialize ordering, and (ii) a unit-level and
  end-to-end reproduction through the real `make_new_reaction` / `check_for_existing_reaction`
  path with isomorphism de-duplication. I did **not** run `rmg.py` on an input that loads a seed +
  overlapping library to observe the discard in a live model + Chemkin file end-to-end; the
  make_new_reaction path is the same code such a run would execute, but the whole-run behaviour
  (including any later step that might re-introduce the sourced reaction) was not observed.
- **The exact library<->family and library<->library collision counts** are run-time quantities and
  were not measured (see Q5).
- **The `--restart` inversion was established from code** (main.py:734 loads the prior core as a
  seed), not exercised in a live restart.
- **Non-neutral (plasma) chemistry** was not the focus: the `get_electron_placement_counts`
  half of the predicate consults placement declarations ungated (including for net-zero
  reactions), which is a separate axis (cf. i134) and could interact with this ordering issue for
  ionised species; that interaction was not probed here.

## Bottom line

"The sourced value wins" is **not a guarantee.** It is true for the library-vs-family case only
because of the fixed registration order, and that same fixed order makes a seed-borne estimate
beat a sourced library value — silently in the common case, and with only a non-committal INFO
"not new" line in the inversion case. The shape is model-wide, not plasma-specific, and has a
second face in library-vs-library collisions.
