# I-102 — quarantining the `Cation_R_Recombination` Marcus rates

Marcus rate rules in `Cation_R_Recombination` evaluate at 10⁻²³ to 10⁻²²⁶ m³/(mol·s) in a
gas-phase or plasma mechanism. They were fitted as electrode electron-transfer rates; the rate law
is being evaluated correctly, but the evaluation has lost the semantics of its source model.
The governing ruling requires them marked, gated with a hard failure, and **preserved**.

This note records what was built, what the evidence shows, and — at the end — what the evidence
could not reach.

---

## 1. The affected set, enumerated from the database

`probes/enumerate_affected.py` loads every family in `database.directory`, asks each whether it
carries a quarantine manifest, and lets the manifest's own criterion select entries out of what is
actually loaded. Nothing in the probe or in RMG knows which entries are affected.

```
families loaded    : 135
families quarantined: 1

family                 : Cation_R_Recombination
state                  : QUARANTINED FOR QUANTITATIVE PLASMA USE
appliesToKineticsClass : Marcus -> <class 'rmgpy.kinetics.arrhenius.Marcus'>
reason                 : electrochemical reference/domain unavailable
rate rules affected    : 7      [1..7]  Root … Root_N-2R->C_N-2BrClFHILiNOPSSi->S
training entries affected: 5    [0..4]  NH2+Li, C2H5+Li, CH3+Li, CH3NH+Li, SH+Li
unaffected in this family: 0 rate rules, 0 training entries
TOTAL AFFECTED ENTRIES : 12
```

**7 rate rules + 5 training entries = 12.** This matches the count in the ticket exactly.
Full output: `repro-evidence/enumerate_affected.log`.

A separate check over the whole database confirms the scope: `Marcus(` appears as a *kinetics
object* in exactly two files, both inside `Cation_R_Recombination`. The other two files that
contain the word "Marcus" — `H_Abstraction/NIST/reactions.py` and `R_Recombination/NIST/reactions.py`
— use it only in prose (`Rice-Ramsperger-Kassel-Marcus`, "the Marcus expression for correlation
between reaction barriers and energetics"), never as a kinetics model. No kinetics **library**
carries Marcus.

---

## 2. The marker: a sidecar manifest that declares a criterion, not a list

`input/kinetics/families/Cation_R_Recombination/quarantine.py` (new file, RMG-database):

```python
name = "Cation_R_Recombination/quarantine"
state = "QUARANTINED FOR QUANTITATIVE PLASMA USE"
appliesToKineticsClass = "Marcus"
reason = "electrochemical reference/domain unavailable"
```

**Why a criterion and not a list.** The ticket names the failure mode directly: "a hard-coded list
in the code that can drift out of step with the database". The manifest never names an entry. The
affected set is *computed* by testing every rate rule and every training entry of the family
against the declared kinetics class, so:

| what someone does to the data | what happens to the quarantine |
|---|---|
| adds an eighth Marcus rate rule | quarantined automatically, no edit |
| refits a rule to `Arrhenius` | released automatically, no edit |
| renumbers, relabels, reorders entries | nothing changes |
| adds a Marcus entry to a *different* family | **not** covered — that family needs its own manifest |

The last row is the honest limit, and §7 below says what it would take to close it.

**Why in the database and not in the code.** The marker belongs with the data it describes.
Quarantining another family is one new file and zero code changes; lifting a quarantine is
deleting one file. There is no registry in RMG-Py enumerating quarantined families — which is why
`test_no_hard_coded_entry_list_exists_in_the_code` can assert that `rmgpy/data/kinetics/quarantine.py`
contains no family name, no node label, and no entry index at all.

**Why `.py` and not `.json`.** Every other file in a kinetics family directory is a Python-literal
data file loaded by `exec` with builtins stripped, and this one is loaded the same way. Three JSON
files exist under `input/`, all ML blobs under `thermo/sidt`; a JSON manifest in a family directory
would be the odd one out for a reviewer.

**Loud on typos.** `appliesToKineticsClass` is resolved against `rmgpy.kinetics` when the family
loads. A misspelling raises `DatabaseError` rather than producing a manifest that announces a
quarantine while gating nothing — which would be worse than no manifest, because it reads as
protection.

---

## 3. Where the gate belongs — and why the ruling's wording moves it

The ruling says hard-fail when a rule "reaches the core mechanism". **Probing that premise first
changed the design.** Running the existing plasma preflight deck against the quarantine database
worktree, before any code change (`repro-evidence/RED-before.stdout.log`):

```
Added 0 new core reactions
Created 1 new edge reactions
    [Lip](3) + [CH3](4) <=> CH3Li(2)
After model enlargement:
    The model core has 6 species and 0 reactions
    The model edge has 2 species and 1 reactions
```

The Marcus reaction lands in the **edge**, never the core. The run does not stop there; it runs on
and dies much later in the Chemkin writer, and only because Chemkin cannot represent Marcus at all
(`repro-evidence/RED-before.stderr.log`). Two consequences:

1. **A gate placed literally at `add_reaction_to_core` would not have fired on this deck.** The
   rate would sit in the edge at ~10⁻²²⁶, which is exactly how a real channel gets silently
   declared unimportant and never promoted. That is a silent outcome of the kind §1 forbids.
2. **The only thing standing between these rates and a quantitative mechanism today is an
   export-format accident.** Any consumer that can represent Marcus integrates them without a
   murmur — RMS/Julia can (`rmgpy/rmg/reactionmechanismsimulator_reactors.py:660`,
   `rmgpy/yaml_rms.py:172`).

So the gate goes at **model admission**, upstream of both edge and core:

- **Primary** — `CoreEdgeReactionModel.apply_kinetics_to_reaction` (`rmgpy/rmg/model.py`). This is
  the single call site where family-estimated kinetics meet a reaction, and the one place where the
  full `(kinetics, source, entry, is_forward)` provenance tuple is in hand. The check runs
  *before* the direction flip and *before* `reaction.kinetics = kinetics`, so a refused reaction is
  left exactly as generated.
- **Backstop** — `add_reaction_to_core` and `add_reaction_to_edge`, for the paths that never pass
  through estimation: seed mechanisms, reaction libraries, and API callers using
  `generate_kinetics=False`.

**How this leaves the electrode/electrolyte door open.** The gate is inside RMG's mechanism
builder and nowhere else, strictly after generation. Deliberately untouched:

| still works | test |
|---|---|
| the database loads | `TestTheDoorLeftOpen::test_the_database_still_loads` |
| the family stays registered and still generates reactions | `…::test_the_family_is_still_registered_and_still_generates_reactions` |
| `Marcus.get_rate_coefficient` still evaluates | `…::test_the_rate_law_still_evaluates` |
| `Reaction.get_rate_coefficient` — not touched at all | pre-existing `reactionMarcusTest.py` still green |
| the rules and training entries are readable and refittable | `TestTheRealQuarantinedFamily::test_nothing_was_deleted` |

A consumer that supplies the missing electrochemical reference gets the same unmodified data.
Only the quantitative RMG mechanism build refuses.

**Cost to ordinary chemistry:** 308 ns per admission for a family with no manifest — about 31 ms
across a 100 000-reaction edge (`repro-evidence/gate-overhead.log`).

---

## 4. The hard failure, with all five required fields

`repro-evidence/AFTER-hardfail.stderr.log`, exit code 1, on the same deck that previously ran on:

```
rmgpy.exceptions.QuarantinedKineticsError: QUARANTINED FOR QUANTITATIVE PLASMA USE: refusing to
admit a quarantined rate at kinetics estimation for the reaction model.
  family:             Cation_R_Recombination
  provenance:         Cation_R_Recombination/training entry 2 "CH3 + Li <=> CH3Li"; rank 3;
                      template [Root_2R->C]; comment: Matched reaction 2 CH3 + Li <=> CH3Li in
                      Cation_R_Recombination/training This reaction matched rate rule
                      [Root_2R->C] family: Cation_R_Recombination
  reaction:           [Lip](3) + [CH3](4) <=> CH3Li(2)
  kinetics class:     Marcus
  reason:             electrochemical reference/domain unavailable
  quarantine record:  …/Cation_R_Recombination/quarantine.py
This reaction is not being dropped, averaged, zeroed, reversed, made irreversible, or given a
substitute rate -- any of those would let the run report success on a mechanism that is quietly
wrong. The family and its data remain in the database untouched. Remove …/quarantine.py only when
the quarantine is genuinely resolved.
```

All five: family ✓, provenance ✓ (training entry 2, rank 3, node `Root_2R->C`), generated reaction
✓, kinetics class ✓, reason ✓.

The run now stops at the **first generated reaction** — before the edge, before any integration,
before any export. Before the change it reached the Chemkin writer.

The provenance field has three sources of decreasing quality and no fourth: the `Entry` when the
estimate was an exact match, the kinetics `comment` when it was averaged (`entry` is `None` there),
and the reaction's template as a last resort. `test_provenance_is_never_blank` pins all three, so
the refusal cannot silently degrade to four fields.

---

## 5. The seven forbidden silent behaviours

`test/rmgpy/data/kinetics/quarantineTest.py::TestTheSevenForbiddenSilentBehaviours` — one test per
bullet, all live-path (they exercise the real `CoreEdgeReactionModel` methods, not a mock of the
gate).

| # | forbidden | how it is tested |
|---|---|---|
| 1 | remove the reaction after generation | the gate raises; it has no sentinel return a caller could read as "skip this one"; the reaction object comes out unmutated and `core`, `edge`, and `new_reaction_list` are all empty |
| 2 | substitute an average rule | RMG's *own* averaged estimate returns `entry=None`, so a gate keyed on exact provenance would let exactly this case through. Both the averaged and the exact-match forms are asserted to refuse identically |
| 3 | evaluate it with `potential = 0` | the refusal precedes binding, so `reaction.kinetics` is still `None` afterwards and there is nothing for anything to evaluate; `check_quarantine` takes no potential parameter, asserted by introspection |
| 4 | mark it irreversible and continue | `reaction.reversible` is `True` before and after; the exception leaves the method (no "continue") |
| 5 | derive a reverse rate | **the live case, not a hypothetical** — this deck reaches the family through its auto-derived reverse template, so `is_forward=False` and RMG's next act is to swap reactants and products. The gate sits before the flip; reactants and products come out in the order they went in |
| 6 | seed a product or cation to bypass it | seeding is precisely how a reaction enters the model *without* kinetics estimation. Both `add_reaction_to_core` and `add_reaction_to_edge` refuse a pre-attached quarantined rate |
| 7 | replace it with a generic collision rate | the gate's only non-raising outcome is `None`; it returns no kinetics at all, and the refused reaction keeps the kinetics it had — none |

**Every one of these was confirmed RED first.** `probes/no_gate_plugin.py` replaces
`rmgpy.rmg.model.check_quarantine` with a no-op, reproducing the pre-ticket code path exactly.
Under it (`repro-evidence/RED-before-tests-no-gate.log`): **10 failed, 26 passed** — all seven
forbidden-behaviour tests and all three gate-placement tests fail. The 26 that still pass are
manifest-loading and door-left-open tests, which are honestly not gate-dependent.

Test 5's RED output is worth reading: with the gate removed, `apply_kinetics_to_reaction` gets far
enough to swap reactants and products (the repr shows `CH3Li` promoted to reactant) before dying
on an unrelated `pairs=None`. The forbidden reverse-rate derivation is not hypothetical; it is what
the code does without the gate.

---

## 6. Negative control — ordinary chemistry is numerically unchanged

`examples/rmg/minimal`, run against the same database before and after
(`repro-evidence/NEGCTL-diff.log`):

```
  BEFORE  The final model core has 26 species and 66 reactions
  BEFORE  The final model edge has 146 species and 332 reactions
  AFTER   The final model core has 26 species and 66 reactions
  AFTER   The final model edge has 146 species and 332 reactions

  IDENTICAL  chemkin/chem.inp
  IDENTICAL  chemkin/chem_annotated.inp
  IDENTICAL  chemkin/chem_edge_annotated.inp
  IDENTICAL  chemkin/species_dictionary.txt

2850089d64e3d74423d22bc26ad2b1969741c93256ebd00beacb34a99cdb54e5  negctl-before/chemkin/chem_annotated.inp
2850089d64e3d74423d22bc26ad2b1969741c93256ebd00beacb34a99cdb54e5  negctl-after/chemkin/chem_annotated.inp
```

Byte-identical, core and edge.

---

## 7. What this could not reach, and what would close it

**The manifest covers a family, not a defect class.** If Marcus data with the same missing
electrochemical reference is ever added to a *different* family or to a kinetics library, it is not
quarantined until someone drops a manifest there. No mechanism currently exists to say "this class
of data, wherever it appears" — the criterion is scoped to the family that owns it. Closing that
would need a database-level manifest (`input/kinetics/quarantine.py`) with the same
criterion-not-list shape, and a loader hook at `KineticsDatabase.load` rather than
`KineticsFamily.load`. That is a real extension, not a refactor, and it is not built here. What
holds the line meanwhile is §1's enumeration: it is a whole-database sweep, so a Marcus entry
appearing elsewhere shows up as an *unquarantined* family the moment anyone runs the probe.

**Kinetics libraries cannot carry a manifest at all.** `load_family_quarantine` is called from
`KineticsFamily.load`; `KineticsLibrary` has no equivalent hook. No library in this database
carries Marcus, so nothing is currently escaping, but the asymmetry is real.

**The electrode/electrolyte consumer was not exercised end to end.** The one in-repo deck that is
plausibly such a consumer, `examples/rmg/SEI_pure_ACN/input.py`, lists `Cation_R_Recombination` and
uses `liquidSurfaceReactor`. It **cannot be run in this environment at all**: it needs
ReactionMechanismSimulator, and `juliacall` is not installed, so it dies in
`reactionmechanismsimulator_reactors.py` while reading the input file, before any family generates
anything — before and after this change alike
(`repro-evidence/SEI-consequence.stderr.log`). So the following is **analysis, not measurement**:
if that deck generates a reaction from `Cation_R_Recombination`, it will now hard-fail, because the
gate is unconditional within an RMG mechanism build.

That is deliberate and follows the ruling, which makes hard-failure categorical and offers no
opt-out — an exemption keyed on reactor type would be a way to make a run appear to succeed, and
the ruling forbids the whole family of those. But it is a real consequence for an electrode model
built *inside RMG*, as distinct from one that consumes the family's data directly, and it is the
one place where "preserve the family for an electrode/electrolyte model" and "the workflow must
hard-fail" genuinely pull against each other. **Flagging it rather than deciding it**: relaxing the
gate for a declared electrochemical domain would be a change to the ruling, not to this code.

---

## Files

**RMG-Py**
- `rmgpy/data/kinetics/quarantine.py` — new; manifest loader, criterion, and the gate
- `rmgpy/exceptions.py` — new `QuarantinedKineticsError`
- `rmgpy/data/kinetics/family.py` — loads the manifest into `KineticsFamily.quarantine`
- `rmgpy/rmg/model.py` — the gate call in `apply_kinetics_to_reaction`, plus core/edge backstops
- `test/rmgpy/data/kinetics/quarantineTest.py` — new; 36 tests
- `docs/i102-quarantine/` — probes, deck, and evidence

**RMG-database**
- `input/kinetics/families/Cation_R_Recombination/quarantine.py` — new; the only change
