# i209 — "the sourced value wins": a guarantee, or an accident of ordering?

Read-only investigation on branch `i209-priority` (`plasma` HEAD `4910d8c52`). No `rmgpy/`
behaviour was changed. Deliverables: this report, `docs/i209-source-priority/repro.py`, and
`test/rmgpy/i209SourcePriorityTest.py` (6 characterization tests, pass on unmodified `plasma`).

This is a revision. The core finding stands and is unchanged: RMG decides "same reaction" from
species references + a per-side electron-placement count only (kinetics and provenance ignored),
keeps whichever reaction was registered first, and does so effectively silently. What changed is
the **headline reachability example**: the previous draft's `--restart`/seed story does **not**
invert the guarantee (Correction 1 below, verified against the code); the durable, general
inversion is **edge pruning** erasing the incumbent (Correction 2 below, confirmed empirically).

## Verdict (Q3), verbatim

**The estimate can win, under these conditions:** a reaction library is loaded, and mid-run RMG's
edge **pruning** removes an edge species that participates in one of that library's reactions.

What the tests confirm — the two steps I directly observed:
1. A prune deletes the sourced library reaction from RMG's global reaction registry
   (`remove_species_from_edge`, def model.py:1515; its deletion block model.py:1578-1591 walks
   *every* `reaction_dict` key, `KineticsLibrary` keys included, and removes entries for the pruned
   species). Confirmed: after `remove_species_from_edge`, the incumbent is no longer in
   `reaction_dict`.
2. A subsequent duplicate check over the same species then finds **no incumbent**
   (`check_for_existing_reaction` returns `False`), so a re-proposed family estimate faces no
   collision to lose to.

Pruning is a default, routinely-triggered mechanism (`tol_keep_in_edge` rate cutoff and the
`maximum_edge_species` cap, model.py:1471-1476), needs no seed mechanism and no special chemistry.
So "the library loads first, therefore the sourced value always wins" is **not durable**: first
registration is erasable, and after erasure nothing protects the sourced value from a regenerated
estimate.

**Named gap — the unobserved join.** I did not run the full sequence end to end in a live model:
prune the incumbent → family enlargement actually re-introduces the same species and regenerates
the same reaction → the estimate is installed as the operative rate in a running mechanism. Steps 1
and 2 are each confirmed directly (see CASE2 / `test_edge_pruning_erases_a_library_incumbent`), but
the joined live prune-then-regenerate-then-install sequence was not exercised; it requires a full
`rmg.py` run (see "What the probe could not reach").

## The mechanism (Q1) — unchanged, confirmed from code and from a run

`are_identical_species_references` (rmgpy/rmg/model.py:2356) decides "same reaction" from **only**
the reactant/product species references (either direction, plus `specific_collider`) and a per-side
electron-placement count from `get_electron_placement_counts`. For every neutral reaction that
count is `(0, 0)` on both sides, so the verdict reduces to a pure heavy-species reference
comparison — kinetics, provenance, and measured-vs-estimated status are never read. The function's
own docstring states the consequence: *"the model then silently keeps whichever was offered
first."* `Species.__eq__` is **reference equality** (`self is other`, species.py:209-211), so the
comparison is by object identity; distinct species only compare equal after `make_new_species` has
unified them (a fact that matters for library-vs-library collisions, Q5).

`make_new_reaction` (model.py:536) keeps the incumbent:
```python
        if check_existing:
            found, rxn = self.check_for_existing_reaction(forward)
            if found:
                return rxn, False          # model.py:599-602 — returns the INCUMBENT, drops `forward`
```

**Which one the model keeps, and why.** Whichever reaction was registered first over those species,
regardless of kinetics or owner. In a normal run the reaction library is registered
(main.py:795-796) before family enlargement (main.py:867-872), so the *library* reaction is the
incumbent and the family estimate is the dropped newcomer — the model keeps the **sourced library
value**. This is decided purely by registration order, not by any preference for a sourced rate.

### Reproduction (exact code and output)

`docs/i209-source-priority/repro.py`, run in `rmg_env` against the committed
`test/rmgpy/test_data/testing_database` (families `H_Abstraction`, `Disproportionation`; libraries
`GRI-Mech3.0`, `ethane-oxidation`). Species `[H] + C=C[CH2]C <=> C=C=CC + [H][H]`; sourced
Arrhenius `A = 1.0e12 cm^3/mol/s` (`1.0e6` SI), estimate `A = 9.9e9 cm^3/mol/s` (`9900` SI). CASE2
and CASE3 go through the real `make_new_species` dedup path (species enter via `make_new_species`;
duplicates resolve through `make_new_reaction` / `check_for_existing_reaction`, over fresh species
objects in CASE3); CASE1 does **not** dedup — it shares the same species objects across both
reactions and calls `register_reaction` / `check_for_existing_reaction` directly, isolating the
predicate's order-dependence. None of the three needs a full `rmg.py` run.

```
CASE1a library-first: are_identical=True found=True kept_A=1000000.0000000001 (SOURCED=1000000.0000000001)
CASE1b estimate-first: found=True kept_A=9900.000000000002 (ESTIMATE=9900.000000000002, SOURCED discarded)
CASE2 before prune: registry_key=GRI-Mech3.0 found=True
CASE2 after prune:  registry_key=None found=False  (incumbent erased; a re-proposed estimate would face no collision)
CASE3 confirmed library-vs-library collisions (GRI first, ethane-oxidation offered): 14 of 18
```

CASE1 shows the order-dependence (same predicate, opposite outcome). CASE2 is the headline: the
sourced library reaction is registered under the `GRI-Mech3.0` key, is the incumbent
(`found=True`), and after pruning species `C=C[CH2]C` it is **gone from the registry**
(`registry_key=None`, `found=False`) — a re-proposed estimate would face no collision. CASE3
confirms library-vs-library collisions through the identity path (Q5).

`test/rmgpy/i209SourcePriorityTest.py` pins all of this (6 tests, all passing on unmodified
`plasma`): `test_edge_pruning_erases_a_library_incumbent` is the headline; the order-dependence,
end-to-end `make_new_reaction`, and identity-confirmed library-vs-library cases are the others.

## Correction 1 — the `--restart`/seed route does NOT invert (I verified; the reviewer is right)

**Independent conclusion: CONFIRMED WRONG — my previous `--restart` headline does not hold.** An
auto-generated seed mechanism does not reload a family estimate as a library reaction that later
shadows a real library. The chain, checked at the exact lines:

- **library.py:356-372** — in `get_library_reactions`, when the library `self.auto_generated` and an
  entry's `long_desc` contains `'rate rule'`, the entry is reconstructed as a **`TemplateReaction`**
  with `family` parsed from the comment — not a `LibraryReaction`. (The earlier branch,
  library.py:343-355, rebuilds "Originally from reaction library" entries as `LibraryReaction`s —
  those are sourced, not estimates.)
- **model.py:1743-1751** — in `add_seed_mechanism_to_core`, a `TemplateReaction` is converted to a
  `LibraryReaction` **only** `if ... not (rxn.family in family_names)`, i.e. only when its family is
  **not** loaded. With the family loaded (the normal case) it stays a `TemplateReaction` and is
  registered under the **family** key.
- **model.py:514-516** — the seed sweep inside `check_for_existing_reaction` iterates
  `for library in self.reaction_dict:` but only descends into keys where
  `isinstance(lib_obj, KineticsLibrary) and library != rxn.family`. A family-keyed
  `TemplateReaction` (its owner is a `KineticsFamily`) is therefore **never** in that shortlist, and
  the first-loop shortlist is keyed under the *newcomer library's* own family, so it is not there
  either. Net: the later real-library reaction is **not** found as a duplicate and is **not**
  discarded.

So run-1 family estimates beating run-2 library rates via `--restart` **does not occur** by that
mechanism. (Refinement, from reading the same code: a **user-supplied, non-auto-generated** seed
mechanism takes the library.py:373 `else` branch and becomes `LibraryReaction`s under a library
key, which the seed sweep *does* consult — so a hand-provided seed whose rates are rough, loaded
alongside a library covering the same reaction, still inverts by ordering. That is a weaker
"estimate vs sourced" story than pruning and is not the headline; the pruning path supersedes it.)

## Correction 2 — edge pruning erases library incumbents (I verified; this is the headline)

**Independent conclusions:**

**Does pruning delete library keys?** **Yes.** `remove_species_from_edge` (def model.py:1515), in
its deletion block model.py:1578-1591, contains:
```python
        for family in self.reaction_dict:                       # model.py:1578 — EVERY key, unfiltered by owner type
            if spec in self.reaction_dict[family]:
                del self.reaction_dict[family][spec]
            for reactant1 in self.reaction_dict[family]:
                if spec in self.reaction_dict[family][reactant1]:
                    del self.reaction_dict[family][reactant1][spec]
            for reactant1 in self.reaction_dict[family]:
                for reactant2 in self.reaction_dict[family][reactant1]:
                    ...
                    for templateReaction in self.reaction_dict[family][reactant1][reactant2]:
                        if spec in templateReaction.reactants or spec in templateReaction.products:
                            ... remove
```
The loop is **not** filtered to `KineticsFamily` keys — it walks every `reaction_dict` key, and
library reactions are registered under their library-name key (`register_reaction` keys on
`rxn.family`, which for a `LibraryReaction` is the library name). The innermost variable is named
`templateReaction` but it removes any reaction object in the bucket, `LibraryReaction`s included.
CASE2 in `repro.py` confirms this empirically: a reaction registered under the `GRI-Mech3.0` key is
present before the prune (`registry_key=GRI-Mech3.0`) and absent after
(`registry_key=None`, `check_for_existing_reaction` returns `False`).

**What triggers the prune?** `CoreEdgeReactionModel.prune` (model.py ~1420-1500). A species is added
to `species_to_prune` when either (i) its `max_edge_species_rate_ratios[index] < tol_keep_in_edge`
(model.py:1471 — the flux cutoff), or (ii) the edge is over the `maximum_edge_species` cap
(model.py:1475-1478, pruning the lowest-rate species while
`num_prunable_species - len(species_to_prune) > maximum_edge_species`, guarded by `tol_move_to_core`).
Both branches then call `remove_species_from_edge` (model.py:1499 and 1507). Both are ordinary,
on-by-default model-enlargement controls.

**Reachability caveat I could not fully close without a run:** pruning removes an *edge* species and
all its reactions, so the inversion completes only if the same chemistry is later regenerated (the
species re-enters via a family acting on core species and the family reaction is re-proposed with no
incumbent). I demonstrated the two load-bearing facts directly (the library incumbent is erased;
a re-proposed reaction then faces no incumbent), but I did **not** run a full `rmg.py` to observe a
prune-then-regenerate sequence end to end in a live model. See "What the probe could not reach".

## Ordering (Q2) — unchanged, file:line that settles it

`rmgpy/rmg/main.py`, `RMG.initialize`: **790-791** seeds -> core; **795-796** libraries -> edge;
**867-872** the `enlarge(spec)` family-generation calls, after both. Registration order is fixed
(**seeds -> libraries -> family reactions**) and the two loops are unconditional and hard-coded; an
input file cannot interleave or reorder them. So the library-beats-family property is genuinely
guaranteed by ordering *at first registration* — but the pruning path (Correction 2) shows that
guarantee is not durable across a run.

## Loudness (Q4) — unchanged: asymmetric, never names the substitution

- **Family reaction discarded (a library/seed was the incumbent) — SILENT.**
  `process_new_reactions` (model.py:972-975): `if not is_new: … continue`. No log at any level.
- **Library/seed reaction discarded — INFO only.** `add_seed_mechanism_to_core` (model.py:1760-1761)
  and `add_reaction_library_to_edge` (model.py:1888-1889):
  `if not isNew: logging.info("This library reaction was not new: {0}".format(rxn))`. INFO-level,
  says only "not new", never that a measured rate was dropped for an estimate, does not name the
  survivor, does not compare kinetics.
- **Pruning itself** logs `logging.info("Pruning species %s", spec)` (model.py:1497, 1505) — it announces
  the species removal, but there is **no** message that a *sourced library reaction* was deleted
  along with it; the erased reaction and its rate simply vanish from the registry. Chemkin output
  and the final model report show no diff: the surviving reaction is written with its own kinetics,
  the discarded one never appears. No warning, no diff.

## Blast radius (Q5)

**Method.** Over the plasma-lineage database `/home/alon/Code/RMG-database-plasma/input/kinetics`
(this worktree's `rmgrc` target), source files only, no generation run. Library reactions: one =
one `^entry(` at column 0 across every `reactions.py` under `libraries/`
(`find libraries -name reactions.py -exec grep -h '^entry(' {} \; | wc -l`). Families: recounted as
`find families -mindepth 1 -maxdepth 1 -type d ! -name __pycache__ | wc -l` and cross-checked with
`ls -d families/*/ | grep -v __pycache__ | wc -l`; every one of those directories contains a
`groups.py`.

**Counts.** 190 `reactions.py` files (71 standalone + 119 nested in 10 container libraries) holding
**35,216 library reaction entries**, against **140 reaction families** (corrected from 138 — see
below). Largest single libraries: CurranPentane (3069), CF2BrCl (2455), CH3Cl (2155), JetSurF2.0
(2087), Chung_solvation_corrections (1784), JetSurF1.0 (1438), NOx2018 (1321). Genuinely
plasma/ion/Li libraries total **207 entries — ~0.6%**; plasma/Li families are ~16 of 140 (~11%). So
the file-measurable ceiling on the blast radius is **35,216 library reactions x 140 families**,
~99.4% of the library population being mainstream hydrocarbon / combustion / NOx / halogen / sulfur
/ surface chemistry. The collision shape is **general, not Lithium/plasma-specific**; it was merely
*noticed* in plasma work.

**Family count correction.** Previous draft said 138; the correct count is **140**. Exact command:
`find /home/alon/Code/RMG-database-plasma/input/kinetics/families -mindepth 1 -maxdepth 1 -type d ! -name '__pycache__' | wc -l`
-> `140` (identical under `ls -d families/*/ | grep -v __pycache__ | wc -l`; all 140 have
`groups.py`).

**Library-vs-library — now CONFIRMED, not a proxy.** `Species.__eq__` is reference equality, so a
genuine collision cannot be read off the files (each library defines its own species objects); it is
resolved only after `make_new_species` unifies them. Driving exactly that path (register
GRI-Mech3.0's reactions first, then offer ethane-oxidation's) yields **14 confirmed identity
collisions of ethane-oxidation's 18 reactions** — e.g. `H + CH4 <=> H2 + CH3`,
`CH3 + CH3 <=> C2H6`, `H + H <=> H2` — every survivor being the first-registered GRI incumbent, its
differently-rated ethane-oxidation copy silently dropped. This is pinned by
`test_library_vs_library_collision_is_identity_confirmed` and reproduced in CASE3. The **mainline
name-proxy figures remain an upper bound, unconfirmed at scale**: 50 of 190 library files contain an
H+O2-type branching reaction, 50 a CO+OH reaction, 90 reference CH3 — these are label/formula
counts, not identity matches, but the phenomenon itself is now confirmed real by the 14 above.

**What this count could NOT reach.** The *actual* number of library<->family and library<->library
collisions in the full database is a run-time quantity — family scope is decided by recipe +
group-tree matching on real molecular graphs during generation, and species identity is established
only while a model is built. 35,216 x 140 is a ceiling on the colliding population, not a count of
pairs that collide.

## What the probe could NOT reach

- **No full RMG generation run was performed.** Evidence is code reading plus unit-level and
  end-to-end reproductions through the real `make_new_species` / `make_new_reaction` /
  `check_for_existing_reaction` / `remove_species_from_edge` paths. In particular I did **not** run
  `rmg.py` to observe a live **prune-then-regenerate** sequence installing an estimate where a
  library value used to be; I confirmed the two load-bearing steps (incumbent erased; re-proposed
  reaction finds no incumbent) directly, but not the whole loop in a running model.
- **The exact library<->family and library<->library collision counts** in the full database are
  run-time quantities and were not measured (see Q5).
- **The `--restart` refutation was established from code** (library.py:356, model.py:1743,
  model.py:514-516), not exercised in a live restart. The user-supplied-seed refinement was reasoned
  from the same code, not run.
- **Non-neutral (plasma) chemistry** was not the focus: the `get_electron_placement_counts` half of
  the predicate consults placement declarations ungated (including net-zero reactions), a separate
  axis (cf. i134) that could interact with pruning for ionised species; not probed here.

## Bottom line

"The sourced value wins" is **not a guarantee.** It holds for library-vs-family only because of the
fixed registration order, and that order also lets a first-registered value be **erased by edge
pruning**: I confirmed directly that a prune deletes a library incumbent from `reaction_dict` and
that a subsequent duplicate check then finds no incumbent — the two steps that make the guarantee
non-durable. The full live sequence that follows (a regenerated family estimate then installed as
the operative rate in a running model) is the plausible consequence but was **not** observed end to
end; it is named as a gap above and in "What the probe could not reach." Whatever trace exists is
only an unrelated "Pruning species …" INFO line. The `--restart` route that the previous draft
headlined does **not** invert (family estimates reload as family-keyed `TemplateReaction`s the seed
sweep skips). The shape is model-wide, not plasma-specific, and has a confirmed second face in
library-vs-library collisions (14 identity-confirmed in a two-library probe).
