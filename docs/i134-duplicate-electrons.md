# I-134 — the duplicate check drops the sink, and it was never about lithium

Evidence labels: **[R]** code (file:line) · **[D]** database (file, entry) · **[M]** measured
(command and output shown) · **[I]** inference, with its basis.

Everything below was run in `/home/alon/Code/RMG-Py-i134-duplicate-electrons`, against
`database.directory = /home/alon/Code/RMG-database-i123b-reaudit/input` — printed at the head of
every probe rather than trusted from configuration [M]. Raw output is in
`docs/i134-duplicate-electrons/evidence/`.

---

## 0. Verdict

The lithium mechanism now carries both channels into the generated model [M]. The repair is three
files and 162 added lines, and it is **not** in the plasma code: the defect was in RMG's notion of
reaction identity, and it reached every charged reaction RMG generates.

Two findings changed the shape of the ticket, and both came from probes that refuted what the
reading of the code suggested:

1. **The blast radius is general, not bounded to the two lithium libraries** (§3). Measured before
   any code changed, and gated on.
2. **Comparing net `electrons` does not fix it** (§4). The obvious repair would have shipped
   looking correct while leaving the lithium pair collapsed.

One thing this repair does **not** reach, measured rather than assumed: a kinetics library
carrying both channels still cannot be loaded (§8). That is pinned as a strict xfail, not papered
over.

---

## 1. The duplicate determination, at source

`CoreEdgeReactionModel.check_for_existing_reaction` [R] `rmgpy/rmg/model.py:424-524` (pre-change
line numbers throughout this section; the repair added one import line above them). It has four
return sites that answer "yes, this reaction is already here", and each is guarded by the same two
predicates:

| site | guard |
|---|---|
| `:471` | `rxn_id == rxn_id0` and `are_identical_species_references` — same-direction, family or library |
| `:485` | `rxn_id == rxn_id0[::-1]` and `are_identical_species_references` — reverse-direction, `KineticsFamily` only |
| `:513` | `(rxn_id == rxn_id0) or (rxn_id == rxn_id0[::-1])` and `are_identical_species_references` — cross-library, forward shortlist |
| `:521` | `are_identical_species_references` alone, no `rxn_id` test at all — cross-library, reverse shortlist |

What those two predicates compare [R]:

- `generate_reaction_id` (`:2292`) — `sorted(get_key(s) for s in reactants)` and the same for
  products, where `get_key(spc)` (`:2338`) is `spc.label`. **Labels only.**
- `are_identical_species_references` (`:2346`) — `rxn1.reactants == rxn2.reactants` and the
  product equivalent, in either direction, plus `specific_collider`. **Object references and the
  collider only.**

`electrons` appears **0 times** in `check_for_existing_reaction`, 0 times in
`generate_reaction_id`, and 0 times in `are_identical_species_references` — measured by
`inspect.getsource` rather than by grep, so it is the running code that is being counted [M]
(`evidence/probe-blast-radius-BEFORE.log`, stages 1 and 2).

RMG's canonical representation keeps the electron out of the participant lists and in the scalar
`Reaction.electrons` [R] `rmgpy/reaction.py:105`. So for two reactions over the same heavy
species, there is nothing left for either predicate to tell apart.

---

## 2. The discard, reproduced on the real mechanism, with the dropped reaction named

`python rmg.py docs/i123-integration/input.py` — the shipped lithium deck, both libraries out of
the union database [M]:

```
Adding reaction library PlasmaElectronImpactIonization to model edge...
    Created 1 new edge reactions
        Li(2) => [Lip](3)

Adding reaction library PlasmaRadiativeRecombination to model edge...
This library reaction was not new: [Lip](3) => Li(2)
    Created 0 new edge reactions

MODEL GENERATION COMPLETED
The final model core has 7 species and 1 reactions
RMG_EXIT=0
```

**The reaction that disappears is the sole entry of `PlasmaRadiativeRecombination`,
`[Lip](3) + e-(1) => Li(2)`** — the cation's only loss channel [M]. Not a count: the log names it,
and the written file confirms which one survived.

`chemkin/chem.inp`, `REACTIONS` block, in full [M]:

```
REACTIONS    KCAL/MOLE   MOLES

Li(2)+e-(1)=>[Lip](3)+e-(1)+e-(1)                   1.490e+18 -0.267    162.150
    TDEP/e-(1)/   ! VoronovEIArrhenius exported as a modified-Arrhenius fit of k(Te) ...

END
```

One reaction. `chem_annotated.inp` attributes it to `PlasmaElectronImpactIonization` and carries no
`PlasmaRadiativeRecombination` line at all [M]. `[Lip](3)` appears in the file only as a product
and in the thermo block — never on a reactant side. The cation is produced and never consumed, the
run exits 0, and the only trace is one `INFO` line that reads like routine deduplication.

**Which return site did it.** Recorded at runtime, by wrapping
`are_identical_species_references` and reading the caller's line number, rather than inferred from
reading the function [M]:

```
    recombination offered next: is_new=False, identity calls at lines [(521, True)]
```

`model.py:521` — the cross-library reverse shortlist, the one guard that does not even test
`rxn_id`. Pure reference comparison, either direction.

---

## 3. Blast radius — GENERAL, measured before the fix

`python docs/i134-duplicate-electrons/probe_blast_radius.py`, full output in
`evidence/probe-blast-radius-BEFORE.log`. The family half uses `Plasma_Electron_Attachment`
applied to triplet dioxygen — a real family, real generation, out of the database, **no lithium in
it** [D][M].

| # | pair | path | collapsed | correct? |
|---|---|---|---|---|
| 1 | ionisation vs radiative recombination | library, `:521` | yes | **NO** |
| 2 | `O2 + e- => O2-` vs its true reverse | family, `:485` | yes | yes — control |
| 3 | `O2 + e- => O2-` vs `O2 => O2- + e-` | family, `:471` | yes | **NO** |
| 4 | `O2 + e- => O2-` vs `O2 => O2-` (neutral) | family, `:471` | yes | **NO** |

Answers to the two questions the brief posed:

**Are family-generated reactions affected?** **Yes.** Rows 3 and 4 are family-generated
`TemplateReaction`s collapsing at `model.py:471`, and row 2 shows `:485` is reached as well. The
blindness is not a property of the cross-library loop; it is a property of the one predicate all
four return sites share.

**Does any forward/reverse pair differing only in electrons collapse?** Yes — but the sharper
finding is row 4: **the collapse needs no electron difference and no plasma library at all.** A
plain neutral reaction whose heavy species coincide with a charged one's is discarded as its
duplicate. Attachment/detachment was the expected second instance; neutral-versus-charged is a
third and much wider one.

This was reported and gated on before any code changed, per the brief's instruction.

### What else consumes the same notion of identity [R]

`are_identical_species_references` has exactly four consumers, all in
`check_for_existing_reaction`. But the same net-only notion lives in `Reaction.is_isomorphic`,
which is consumed by:

- `ReactionModel.merge` (`model.py:129`, `either_direction=True`) — `mergeModels.py` would drop
  the sink the same way;
- `KineticsLibrary.check_for_duplicates` (`library.py:368`) and `convert_duplicates_to_multi`
  (`:400`) — see §8;
- `rmgpy/rmg/pdep.py:927,939`.

---

## 4. Why the obvious fix is wrong — the measurement the repair turns on

The natural repair is to compare `Reaction.electrons`. **It does not work**, and it fails
silently.

`Reaction.is_isomorphic` *already* compared the net scalar before this change [R]
`rmgpy/reaction.py:744` (`self.electrons == other.electrons`) and `:767`
(`self.electrons == -other.electrons`). Run it on the real pair [M]:

```
is_isomorphic(ionisation vs recombination) = True
```

RMG's own canonical identity predicate says they are the same reaction. The arithmetic: the
ionisation carries `electrons = +1`, the recombination `-1`; the reverse test asks
`+1 == -(-1)`, which holds. **Equal and opposite is exactly the relation a genuine reverse pair
has**, and the net scalar cannot distinguish that from this.

They are not reverses [M]:

| | placement `(reactant, product)` | net | order |
|---|---|---|---|
| electron-impact ionisation | `(1, 2)` | `+1` | 2 |
| *its* reverse — three-body recombination | `(2, 1)` | `-1` | 3 |
| radiative recombination | `(1, 0)` | `-1` | 2 |

`(1, 2)` reversed is `(2, 1)`, not `(1, 0)`. Different molecularity, different rate law, different
temperature dependence, a photon in one and a third body in the other.

The pair that separates them already exists in the codebase:
`FAMILY_ELECTRON_PLACEMENT` [R] `rmgpy/electron_placement.py:172`, read by the plasma reactor
through `resolve_electron_placement` and by both export writers through `expand_electrons`. I-113
widened it to two-sided declarations precisely because "a net count is ONE number and placement
needs TWO" [R] `rmgpy/electron_balance.py`. **That same argument turns out to govern identity.**
It had simply never been applied there.

---

## 5. The chosen layer, and why the other two are worse here

The brief offered three. **Layer 1 — the identity comparison includes electron stoichiometry —
chosen**, with the correction that "stoichiometry" must mean per-side placement, not the net
scalar (§4).

Argued from the code's own model of identity rather than from needing both reactions: RMG decides
two reactions are the same when they have the same participants on each side plus the same
collider. The electron *is* a participant; it is stored as metadata for representational reasons,
not because it is not one. Every other participant is compared per side. Comparing the electron per
side makes identity consistent with itself, and it is what makes §3 row 2 keep collapsing while
rows 1, 3 and 4 stop.

**Layer 2 — make the reverse-match test reject a match requiring an electron imbalance — rejected
on measurement, not taste.** It addresses `:485` and `:521`. Rows 3 and 4 of §3 collapse at
`:471`, the *same-direction* branch, which layer 2 leaves untouched. It also leaves
`Reaction.is_isomorphic` — and therefore `ReactionModel.merge` and the library duplicate check —
wrong. It would have made the deck green while the general defect stayed.

**Layer 3 — keep the two distinguishable upstream so the question never arises — rejected because
the boundary does not exist and inventing it is the expensive move.** Upstream of the model builder
the reactions are in canonical form by design: the electron is metadata, and that is the
representation the reactor, both writers and the whole database layer agree on. Making the model
builder see explicit-electron reactions would fork the representation at exactly the seam I-126
spent a ticket un-forking. The one genuine upstream fact — the placement declaration — is what
layer 1 reads.

---

## 6. The repair

Three files, one new function, two consumers of it.

**`rmgpy/electron_balance.py`** — `get_electron_placement_counts(reaction)`. Returns
`(reactant_count, product_count)` for the reaction in its own stored orientation. Reads the owner's
declaration where there is one, deriving the orientation from the sign of `electrons` against the
declared net; falls back to the same net-derived rule `expand_electrons` uses where there is not.
The docstring carries §4's argument in full, because the next person to look at this will reach for
`reaction.electrons` and needs to know why it is not equivalent.

**`rmgpy/rmg/model.py`** — `are_identical_species_references` compares those counts side for side
on a same-direction match, and crosswise on a reverse match. All four return sites are repaired by
the one change, because all four route through this predicate.

**`rmgpy/reaction.py`** — `Reaction.is_isomorphic` does the same, replacing its two net-scalar
tests. Cythonized; `make build` run, and the compiled `.so` confirmed to carry the helper [M].

No species, library, family or reaction is named anywhere in the change. The duplicate check is not
weakened: §3 row 2 and §7's duplicate tests are the evidence.

### What can and cannot move, stated exactly

The two predicates were **not** in the same starting state, and conflating them would overclaim:

- `Reaction.is_isomorphic` already compared the net scalar. For an owner with no declaration, the
  per-side counts *reduce* to precisely those two comparisons — `(r1,p1)==(r2,p2)` iff
  `e1==e2`, and `(r1,p1)==(p2,r2)` iff `e1==-e2`. Its verdicts are therefore unchanged **bit for
  bit**, and `TestTheUndeclaredOwnerContract` pins that over the whole grid of electron counts
  from −3 to +3 [M].
- `are_identical_species_references` compared **no** electron information at all. For a charged
  reaction its verdict genuinely changes — that is the repair. What bounds it: a neutral reaction's
  counts are `(0, 0)`, so any two neutral reactions still compare equal on electrons. All of RMG
  outside the charged families and the plasma libraries is neutral, and is untouched.

A one-sided declaration `(n, 0)` or `(0, n)` places exactly where the net rule would, so those
owners move nothing either. `PlasmaElectronImpactIonization`'s `(1, 2)` is the only two-sided
declaration shipped, and therefore the only owner whose answers can change at all [M].

---

## 7. Verification

**The deck, after** [M] — `python rmg.py docs/i123-integration/input.py`, `RMG_EXIT=0`:

```
The final model core has 7 species and 2 reactions

REACTIONS    KCAL/MOLE   MOLES

Li(2)+e-(1)=>[Lip](3)+e-(1)+e-(1)                   1.490e+18 -0.267    162.150
    TDEP/e-(1)/   ! VoronovEIArrhenius exported as ...

[Lip](3)+e-(1)=>Li(2)                               1.734e+14 -0.801    0.010
    TDEP/e-(1)/   ! BadnellRRArrhenius exported as ...

END
```

Both channels. `chem_annotated.inp` attributes them to the two libraries respectively [M]. The
cation now appears on a reactant side: it has a loss channel. The `was not new` line is gone.

**The probe, after** — `docs/i134-duplicate-electrons/probe_blast_radius.py` exits **0**, all
checks passing, with §3's table inverted except row 2 [M] (`evidence/probe-blast-radius-AFTER.log`).
The probe now asserts the repaired behaviour and stands as a regression artifact; the BEFORE log is
kept beside it as the record of the defect.

**New tests** — `test/rmgpy/i134DuplicateElectronsTest.py`, 179 passed + 1 xfail [M].

- RED before the repair: with the helper present but both consumers reverted to `HEAD` and
  rebuilt, **48 failed, 131 passed, 1 xfailed** [M] (`evidence/new-tests-RED-before-fix.log`).
  Isolating the helper this way — it is additive and nothing calls it — is what makes the red
  attributable to the two consumers rather than to a missing import. The six non-parametrised
  failures are the headline ones:

  ```
  TestChargedFamilyReactionsAreNotConfused::test_the_same_direction_charge_reversal_survives
  TestChargedFamilyReactionsAreNotConfused::test_the_neutral_transformation_over_the_same_heavy_species_survives
  TestChargedFamilyReactionsAreNotConfused::test_the_model_keeps_both_through_the_real_path
  TestLithiumChargeNetworkReachesTheModel::test_both_channels_enter_the_model[True]
  TestLithiumChargeNetworkReachesTheModel::test_both_channels_enter_the_model[False]
  TestLithiumChargeNetworkReachesTheModel::test_the_cation_has_a_loss_channel_in_the_model
  ```

  The remaining 42 are the `test_the_model_predicate_is_tightened_but_only_for_charged_reactions`
  grid, which is red for every charged pair and green for the neutral corner — the shape of §6's
  claim, measured.
- The sink is asserted **present in the model** — `CoreEdgeReactionModel.make_new_reaction`
  returning `is_new=True`, in **both** offer orders, since the drop was order-dependent — and
  separately as chemistry: every cation produced is also consumed.
- The non-lithium family pair is pinned in `TestChargedFamilyReactionsAreNotConfused`.
- The duplicate check is pinned intact in `TestGenuineDuplicatesStillCollapse`: identical neutral
  reactions, a neutral reaction and its reverse, identical charged reactions, a charged reaction
  and its true reverse, and collider mismatch.

**Full unit suite, serial** (`python -m pytest -m "not functional and not database"`, no `-n`) [M]:

| | failed | passed | skipped | deselected |
|---|---|---|---|---|
| before | 2 | 2834 | 50 | 131 |
| after | **2** | **3002** | 50 | 143 |

The same two failures before and after, both **pre-existing on this base and unrelated to this
change** [M]:

- `electronPlacementTest.py::test_declaration_registry_is_explicit_and_closed` — the registry
  assertion was not updated when `PlasmaRadiativeRecombination` was added to
  `FAMILY_ELECTRON_PLACEMENT`; it asserts a three-entry dict against the shipped four.
- `preflightDeckFamilyExclusionTest.py::...[docs/i102-quarantine/input.py]` — the quarantine deck
  declares the quarantined family, which is that deck's purpose.

Neither is mine to fix here, and neither moved. **No non-plasma behaviour changed**, so the
brief's stop-and-report condition on that did not trigger.

**Negative control** [M] — `examples/rmg/minimal`, run before and after:

```
before: core 26 species / 66 reactions, edge 146 / 332
after:  core 26 species / 66 reactions, edge 146 / 332
```

`diff -r` over every generated artifact (excluding `*.log`, `*.png`, `*.xls`, `*.h5`): every
Chemkin file, both species dictionaries, `tran.dat`, all 22 RMS files, the seed mechanism, and all
22 solver CSVs are **byte-identical**. Two files differ, `cantera_from_ck/chem.yaml` and
`chem_annotated.yaml`, and the entire diff of each is one line:

```
< date: Thu, 27 Aug 2026 00:10:19 +0300
> date: Thu, 27 Aug 2026 00:26:16 +0300
```

The conversion timestamp. Identical model, identical output files.

---

## 8. Adjacent landmine — measured, and **not** fixed by this change

`KineticsLibrary.check_for_duplicates` runs inside `load` [R] `library.py:599` and raises
`DatabaseError` on any unmarked duplicate, deciding duplication with `Reaction.is_isomorphic`. It
looked as though repairing that predicate would let one library carry both channels. **It does
not** [M]:

```
RESULT: check_for_duplicates REFUSED -> DatabaseError
   Unexpected duplicate reaction [Li] => [Lip] in kinetics library PlasmaBothChannels.
   Reaction index 2 matches index 1.
```

The reason is upstream of anything this change touches. The check runs over `entry.item`, which
`KineticsLibrary.load_entry` builds as a plain `Reaction` **carrying no owner at all** [M]:

```
PlasmaElectronImpactIonization  entry.item  family=None  electrons=+1  decl=None  counts=(0, 1)
PlasmaRadiativeRecombination    entry.item  family=None  electrons=-1  decl=None  counts=(1, 0)
```

No owner means no declaration means the net rule — under which the two are `(0, 1)` and `(1, 0)`,
exact mirrors, and still duplicates. The owner only appears later, in `get_library_reactions`,
which wraps each entry in a `LibraryReaction` whose `family` is the library label; from there on
the declaration is visible and this repair applies [M].

The verdict is also **unchanged** by this commit rather than broken by it: the old code compared
`+1 == -(-1)` (true), the new code compares `(0,1)` against the mirror of `(1,0)` (equal, true).
Same answer, same reason.

Fixing it needs a design decision this ticket has no mandate for: the placement registry is keyed
one declaration per owner, so a single library cannot declare two different placements even once
the owner is recorded. Pinned as `@pytest.mark.xfail(strict=True)` —
`test_one_library_carrying_both_channels_can_be_loaded` — so it converts to a failing test the
moment someone lands it.

**Practical consequence for this campaign:** the two channels must stay in two libraries. The deck
already does this, and it works.

---

## 9. What this pass could not reach

- **`ReactionModel.merge` (`model.py:129`) is repaired but not exercised.** It consumes
  `is_isomorphic(either_direction=True)`, so the fix reaches it by construction, but I ran no
  `mergeModels.py` over two plasma models. **UNKNOWN** whether merge has other electron-blind steps
  of its own. Running `scripts/mergeModels.py` on two decks each carrying one channel would settle
  it.
- **`rmgpy/rmg/pdep.py:927,939`** likewise consume the repaired predicate, untested here. Plasma
  reactions are not pressure-dependent in any shipped deck, so this is latent. **UNKNOWN.**
- **The `check_template_rxn_products` shortcut** [R] `rmgpy/reaction.py:708-714` still compares the
  net scalar, negating on `is_forward`. It is a products-only generation shortcut and I did not
  change it. Whether its `is_forward` negation is even correct is doubtful — `_create_reaction`
  stores *both* directions in family-forward orientation and does not negate — but that is a
  separate question and touching it was outside this repair. **UNKNOWN.**
- **Performance.** The added comparison is O(1) — a `getattr`, a dict lookup, and a lazy import
  statement — on RMG's hottest path. The negative control's own reported execution time was 28 s
  before and 30 s after, and the unit suite 258 s versus 266 s — but the after-suite runs 168 more
  tests, and neither figure is a controlled benchmark on a quiet machine. **No material slowdown
  observed; not measured properly.**
- **Functional and database-marked suites** were not run in full; only the database-marked tests in
  the new file were. The brief asked for the unit suite.
- **The Chemkin read-back defect** (`docs/i123b-reaudit.md` §7.3) is untouched by design — the
  `TDEP/` line still blocks `ck2yaml`, and the deck's Cantera translation step still fails while
  `rmg.py` exits 0. A separate ticket owns it.
- **Regression-test artifacts.** No `test/regression/` entry was added; the new coverage is unit
  and database-marked.

---

## 10. Reproducing this

```bash
cd /home/alon/Code/RMG-Py-i134-duplicate-electrons
conda activate rmg_env
export PYTHONPATH=/home/alon/Code/RMG-Py-i134-duplicate-electrons:$PYTHONPATH
python -c "import rmgpy; print(rmgpy.__file__)"   # must be under this worktree
make build                                         # never bare `make`

python docs/i134-duplicate-electrons/probe_blast_radius.py      # exit 0
python -m pytest test/rmgpy/i134DuplicateElectronsTest.py       # 179 passed, 1 xfailed
python -m pytest -m "not functional and not database"           # 2 failed (pre-existing), 3002 passed

mkdir -p /tmp/i134-deck && cp docs/i123-integration/input.py /tmp/i134-deck/
python rmg.py /tmp/i134-deck/input.py \
  > >(tee /tmp/i134-deck/stdout.log) 2> >(tee /tmp/i134-deck/stderr.log >&2)
sed -n '/^REACTIONS/,/^END/p' /tmp/i134-deck/chemkin/chem.inp   # both channels
```
