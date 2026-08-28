# I-126 — which layer owns the canonical → explicit electron conversion?

**Verdict: layer 1, narrowed to its converter.** The export path resolves, but not in the
writers — in `rmgpy/electron_balance.expand_electrons`, the single function all export writers
already share, which now reads the *same* `FAMILY_ELECTRON_PLACEMENT` declaration the plasma
reactor reads. Layer 2 is **not implementable** without changing `Reaction.is_balanced`, measured
below. Layer 3's boundary **does not exist**: at the only point both the reactor and the writers
are downstream of, the electron species is not yet in the model.

**The guard caught a real defect.** `check_electron_reactant_order` is untouched, and this ticket
argues it was right rather than in the way.

| | |
|---|---|
| Worktree | `/home/alon/Code/RMG-Py-i126-chemkin-electrons`, branch `i126-chemkin-electrons`, base `3559bf431` |
| Database | `/home/alon/Code/RMG-database-i123-integration/input` (`c7bd96292`), pinned by the worktree's `rmgrc`, printed at the head of every run below |
| Files changed | `rmgpy/electron_balance.py` (the fix), `rmgpy/electron_placement.py` (a docstring this ticket refutes), `docs/i123-integration/input.py` (writer config), one new test file |
| Deck | `python rmg.py docs/i123-integration/input.py` — **RED before, exit 0 after** |
| New tests | `test/rmgpy/i126ChemkinElectronOrderTest.py`, 17 passed |
| Negative control | `examples/rmg/minimal/input.py` — `chem.inp`, `chem_annotated.inp`, `chem_edge.inp`, `species_dictionary.txt` all **byte-identical** before and after |

Evidence labels: **[R]** code (file:line) · **[D]** database (file, entry) · **[M]** measured
(command and output shown) · **[I]** inference with its basis · **UNKNOWN** where nothing was
established.

---

## 1. The RED, before anything was changed

**[M]** Base commit `3559bf431`, `make build` exit 0, run from the worktree root so the pinned
`rmgrc` resolves:

```
$ python -c "import rmgpy; print(rmgpy.__file__)"
/home/alon/Code/RMG-Py-i126-chemkin-electrons/rmgpy/__init__.py
$ python -c "from rmgpy import settings; print(settings['database.directory'])"
/home/alon/Code/RMG-database-i123-integration/input

$ python rmg.py -o <scratch>/red-run docs/i123-integration/input.py
...
    The model core has 6 species and 0 reactions
    The model edge has 1 species and 1 reactions
Completed initial enlarge edge step.
Saving current model core to Chemkin file...
Chemkin file contains 0 reactions.
Saving current model core and edge to Chemkin file...
Traceback (most recent call last):
  File ".../rmgpy/rmg/main.py", line 961, in execute
    self.save_everything()
  File "rmgpy/chemkin.pyx", line 1938, in rmgpy.chemkin.write_kinetics_entry
    reaction_string = write_reaction_string(reaction, java_library, species_list=species_list)
  File "rmgpy/chemkin.pyx", line 1829, in rmgpy.chemkin.write_reaction_string
    check_electron_reactant_order(reaction, reactants, reaction_string)
rmgpy.exceptions.MechanismWriterError: Reaction Li(2) => [Lip](3) has kinetics
VoronovEIArrhenius whose rate is proportional to the electron density, but the exported
equation "Li(2)=>[Lip](3)+e-(1)" has no electron among its reactants, so a solver would
evaluate it at the wrong reaction order.
RED_EXIT=1
```

Note where it dies: the **core** save (0 reactions) succeeds and the **edge** save fails. The
ionisation reaction is on the edge at that point, and the plasma reactor has not run yet. So this
is purely an export-path failure — which is why the audit's probe, which never exports, is green
on the same tree.

---

## 2. Did the guard catch a real defect, or is it a false positive in itself?

**It caught a real defect, and the defect is not in the guard.** Four independent lines of
evidence, in ascending order of force.

**2.1 The rate law's own class docstring says so.** **[R]**
`rmgpy/kinetics/arrhenius.pyx:1505-1512`, in `VoronovEIArrhenius`:

> Writing the incident electron out is not decoration. It is what makes the exported equation
> second order, which is what the rate coefficient's `cm^3/(molecule*s)` dimensionality says it
> must be; the metadata electron lands on the product side and cannot supply it. `H => H+` with
> `electrons = +1` balances and loads, but `rmgpy.electron_balance.check_electron_reactant_order`
> refuses to export it, **correctly** — a solver would evaluate it first order in H and never
> multiply by the electron density.

The guard's refusal is a documented, intended property of the rate law, written by whoever
shipped it. It is not an accident of the export path.

**2.2 The two guards disagree with each other about the same file, which is the tell.** **[M]**
Run both export checks on the net-derived expansion of the real library reaction:

```
converter                    equation                            balance      order
expand_electrons             [Li](2)=>[Lip](3)+e-(1)             E-balanced   ORDER REFUSED
resolve_electron_placement   [Li](2)+e-(1)=>[Lip](3)+e-(1)+e-(1) E-balanced   order OK
```

The net-derived equation **balances in `E`**. That is exactly the failure mode
`check_electron_reactant_order` exists for and `check_electron_balance` cannot see: a file that
is charge-correct, element-correct, and dimensionally wrong. Remove the order guard and the
mechanism exports silently.

**2.3 The number.** **[M]** The Chemkin line the writer produces carries
`A = 1.490e+18 cm^3/(mol*s)` — a *second-order* coefficient, because
`_plasma_arrhenius_for_chemkin` reduces the Voronov fit and
`get_conversion_factor_from_si_to_cm_mol_s()` returns `1e6` for it, not `1e0`. Beside a
one-reactant equation, a solver reads that number as `s^-1` and computes `k·[Li]` instead of
`k·[Li]·[e]`. At the deck's own electron density:

```
deck n_e   = 1e17 m^-3 = 1.66054e-07 mol/m^3
1 / [e]    = 6.022e+06 m^3/mol
```

so the pre-guard file would have overstated the ionisation rate by a factor of **6.0 × 10⁶**, at
that density, with nothing in the file to indicate it. That is the "off by a factor of the
electron density" the guard's message names, made concrete.

**2.4 What a correct file looks like without the guard: there isn't one.** The question "show
what a correct exported file looks like without it" has no answer that keeps the canonical
equation. `Li => Li+ + e-` cannot be made correct by any choice of `A`, `n`, `Ea`, because the
rate's dependence on `[e]` is a *functional* dependence, not a constant: rescaling `A` by a fixed
`1/[e]` is correct at one electron density and wrong at every other, and the electron density is
a state variable the mechanism is supposed to solve for. Chemkin's `TDEP` keyword moves the
*temperature* the rate is evaluated at; it does not add a concentration factor. So the guard is
not blocking a correct file behind a formality — there is no correct file on that side of it.

**Conclusion.** The guard is a true positive. It is untouched by this ticket: not relaxed, not
narrowed, not bypassed. What was wrong is the converter that fed it.

---

## 3. The diagnosis: one representation, two converters

RMG carries the electron of a charged reaction as the scalar `Reaction.electrons`, and converts
to the explicit form on the way to each consumer. There are **two** converters, and they were
written against different answers.

| | reactor | export path |
|---|---|---|
| function | `resolve_electron_placement` **[R]** `rmgpy/electron_placement.py:276` | `expand_electrons` **[R]** `rmgpy/electron_balance.py:72` at base `3559bf431`, `:134` now |
| order source | the two-sided `FAMILY_ELECTRON_PLACEMENT` declaration | the net scalar |
| can express `(1, 2)`? | yes | **no** |
| production call sites | `rmgpy/solver/plasma.pyx:236-237, 334` | `chemkin.pyx:1785`, `chemkin.pyx:1935`, `yaml_cantera2.py:917`, `reaction.py:312` |

**[M]** On the real library reaction, loaded from the integration database:

```
canonical reaction : [Li](2) => [Lip](3)   family='PlasmaElectronImpactIonization'
                     electrons=1   kinetics=VoronovEIArrhenius   declaration=(1, 2)

expand_electrons             -> reactants ['[Li](2)']           order 1
                                products  ['[Lip](3)', 'e-(1)']
resolve_electron_placement   -> reactants ['[Li](2)', 'e-(1)']  order 2
                                products  ['[Lip](3)', 'e-(1)', 'e-(1)']
```

**The net scalar under-determines placement, and that is the whole defect.** A net count is ONE
number; placement needs TWO. `electrons = +1` is equally consistent with `Li + e- => Li+ + 2 e-`
(order 2) and `Li => Li+ + e-` (order 1), and `expand_electrons` always produces the second. The
module that owns the declaration says exactly this, at length, in its own docstring **[R]**
`rmgpy/electron_placement.py:58-71` — and then the export path was left reading the one number
that cannot say it.

**Why it only bites now.** Three of the four declared owners are one-sided — `(1, 0)` — and for a
one-sided declaration the net rule and the declaration place electrons **identically**. That
equality is already pinned by test **[R]**
`test/rmgpy/i113IonisationPlacementTest.py:343-365`. `PlasmaElectronImpactIonization` is the only
declaration whose two numbers differ, so it is the only reaction shape on which the two
converters can disagree at all — and the moment the RMG-database ionisation branch put such a
reaction in a model, they did.

---

## 4. The three candidate layers

### Layer 2 — the library entries carry explicit electrons: **NOT IMPLEMENTABLE** as scoped

This is the layer that would remove the two-representations problem at the source, and it is the
one the brief warns to check for a prior decision before overturning. There is one — **[D]**
`input/kinetics/libraries/PlasmaElectronImpactIonization/reactions.py`, `longDesc`: *"the entry
below is the canonical database form and must not be written with an explicit electron
participant (the resolver refuses double representation)."* But the prior decision is not what
blocks it. **The loader is.**

**[M]** Build the layer-2 form in memory and run it through the loader's own acceptance test:

```
hypothetical explicit-electron entry: [Li](2) + e-(1) => [Lip](3) + e-(1) + e-(1)
  is_balanced() with electrons=0    : False
     (Reaction.is_balanced counts the 'e' PSEUDO-ELEMENT by ATOM COUNT:
      left e-atoms = 1, right e-atoms = 2)
  writers' E-count rule instead     : left E=1 right E=1
```

Two separate rules for "does the electron balance", and only one of them can express ionisation:

- **[R]** `rmgpy/reaction.py:1721-1723` — `is_balanced` compares `reactant_elements[element]` to
  `product_elements[element]` for every element including `e`, and returns `False` **before** any
  charge folding. An ionisation with a surplus electron has 1 on the left and 2 on the right, so
  it can never balance under this rule, whatever `electrons` is set to.
- **[R]** `rmgpy/electron_balance.py:228-242` — `get_species_electron_count`, the rule both
  *writers* use, counts `E` as minus the net charge. Under that rule the same equation balances
  (1 = 1), which is why the reactor's view passes step 10 of the resolver.

**[R]** `rmgpy/data/kinetics/library.py:585` raises `DatabaseError: Reaction ... was not
balanced! Please reformulate.` on the first rule. So the entry would not load.

And it cannot be dodged by setting `electrons = 0` in the entry: **[R]**
`rmgpy/data/kinetics/library.py:567-583` copies `entry.data.electrons` onto the reaction
*unconditionally* whenever the rate law has the attribute, and `VoronovEIArrhenius.electrons`
defaults to `+1` and is set in `__init__` ahead of every branch **[R]**
`rmgpy/kinetics/arrhenius.pyx:1535-1537`. **[M]** With that propagation applied, the resolver
then refuses the same entry for a *second* reason:

```
  after rxn.electrons = kinetics.electrons.value = 1:
  is_balanced()  : False
  resolver       : ElectronPlacementError: ... represents the electron twice: 3 explicit
                   electron participant(s) AND a nonzero metadata electron count
```

**Blast radius: three changes across two repositories, one of them a named non-goal.**
`Reaction.is_balanced` (used by every family, every library, every depository, and the reverse-
reaction machinery), `KineticsLibrary.load`'s propagation, and then the database entries. Changing
`is_balanced`'s element rule is not a local edit — it is a change to what "balanced" means for
every charged reaction in RMG, and the ticket names it as a non-goal precisely because it is
load-bearing far outside this failure. **Rejected on measurement, not on taste.**

### Layer 3 — the reaction resolves itself once, at a shared boundary: **NO SUCH BOUNDARY**

The brief asks to establish whether such a boundary exists rather than assume it. It does not,
and the reason is an ordering fact, not a taste.

The reactor and the writers both read `rmg.reaction_model.core.reactions` / `.edge.reactions`.
The only boundary upstream of both is therefore the point a reaction enters the reaction model —
**[R]** `rmgpy/rmg/model.py:1829` `add_reaction_library_to_edge`, `:1590`
`add_reaction_to_core`, `:1623` `add_reaction_to_edge`.

**[M]** From the RED run's own log, in order:

```
Adding reaction library PlasmaElectronImpactIonization to model edge...
    Created 1 new edge species  [Lip](3)
    Created 1 new edge reactions  Li(2) => [Lip](3)
...
Adding species e-(1) to model core          <-- FIVE enlargements later
```

The library reaction enters the model **before the electron species does**. `resolve_electron_
placement` requires exactly one electron species resolvable from the species list and raises
`ElectronPlacementError: No electron species is resolvable from the supplied species list` when
there is none **[R]** `rmgpy/electron_placement.py:499-506`. So resolution at model entry would
hard-fail on the very reaction it exists to serve.

That alone settles it, and there is a second, structural objection worth recording because it
would survive any re-ordering: the model's reaction objects are its *identity*. `PlasmaReactor`
goes out of its way not to leak its resolved views into that identity — it keeps the canonical
lists and re-keys `reaction_index` back to them at the end of `initialize_model` **[R]**
`rmgpy/solver/plasma.pyx:290-303`, with a comment naming the caller (`remove_species`) that
depends on it. Resolving at model entry would put the electron into `Reaction.reactants` for
every consumer of the model: duplicate detection, `make_new_reaction`, flux pairs, pruning,
`react_edge`'s reactant enumeration, and the seed-mechanism writer. **Blast radius: the widest of
the three, for a boundary that cannot run where it would have to run.**

### Layer 1 — the export path resolves: **CHOSEN**, narrowed to the converter

The brief's caution about layer 1 is that "it puts representation logic in a writer, and there is
more than one writer". That caution is answered by *not putting it in a writer*. All four export
call sites already funnel through one function, and that function already lives in the module
whose own docstring says **[R]** `rmgpy/electron_balance.py:41-42`: *"Both writers must behave
identically here, so the logic lives in one place rather than being written twice."*

So the change is to `expand_electrons` itself:

> if the reaction's owner carries a placement declaration, place by the declaration; otherwise
> keep the net-derived rule exactly as it was.

Every export consumer is repaired by one edit, and none of them had to be touched:

| call site | what it needed |
|---|---|
| `chemkin.pyx:1785` | the equation string |
| `chemkin.pyx:1935` | `num_reactants`, which sizes the A-factor unit assertion |
| `yaml_cantera2.py:917` | the equation string |
| `reaction.py:312` (`Reaction.to_cantera`) | the in-memory Cantera object's stoichiometry |

**Blast radius: bounded by the registry, and provably so.** The declared set is four owners
**[R]** `rmgpy/electron_placement.py:172-177`. Three are one-sided, and for a one-sided
declaration the two rules place identically — so the *only* reactions whose exported equation can
change are those of an owner whose declaration has electrons on both sides, which today is
`PlasmaElectronImpactIonization` alone. Everything else in RMG — every charge-transfer reaction,
every reaction with no family attribution, every neutral reaction — takes the same code path it
took before. This is pinned by `TestOneSidedOwnersDoNotMove` and
`TestUndeclaredReactionsKeepTheNetRule`, both of which **passed before the fix as well as after**,
which is what makes them controls rather than decoration.

**The prior argument against this, and why it was wrong.** `_place_declared_electrons`'s docstring
argued that `expand_electrons` must stay purely net-derived because its export callers have "no
family declaration in hand to widen it with", and that widening it would mean "threading a
declaration through four export call sites". **[R]** That premise is false: the declaration is
keyed on `Reaction.family`, which every reaction already carries and which the resolver itself
reads the same way **[R]** `rmgpy/electron_placement.py:323`. Nothing had to be threaded
anywhere. The second half of that argument — *"the export writers' whole purpose is to be the one
place the net scalar IS authoritative"* — is what this failure refutes: for a two-sided shape the
net scalar is not authoritative, it is insufficient, and the guard is what proved it. That
docstring has been corrected in place rather than left contradicting the code.

**What is deliberately *not* done: routing the export path through `resolve_electron_placement`
itself.** It is tempting — same function, guaranteed same answer — and it is wrong here, because
the resolver validates for a *reactor*, not for a file. It refuses pressure-dependent kinetics,
and it refuses any rate law carrying its own `electrons` field outside a two-name allowlist
**[R]** `rmgpy/electron_placement.py:477-484`. `Cation_R_Recombination` is a declared owner whose
reactions can carry `ArrheniusChargeTransfer` kinetics, which has exactly such a field — so
routing export through the resolver would newly refuse to *write* reactions that export correctly
today, in a family this ticket is told not to touch. The declaration is shared; the reactor's
acceptance policy is not. The export path keeps its own guards, `check_electron_balance` and
`check_electron_reactant_order`, which run immediately after and are the file's correctness
criteria.

---

## 5. The other writers — do Cantera and RMS fail the same way?

**Cantera YAML v2: yes, identically, and it is fixed by the same change.** **[R]**
`rmgpy/yaml_cantera2.py:917` calls the same `expand_electrons` and `:935-936` runs the same two
guards. **[M]** Before the fix, on the same reaction:

```
-- Cantera YAML v2 (rmgpy/yaml_cantera2.py) --
   MechanismWriterError: Reaction [Li](2) => [Lip](3) has kinetics VoronovEIArrhenius whose
   rate is proportional to the electron density, but the exported equation
   "[Li](2) => [Lip](3) + e-(1)" has no electron among its reactants
```

Same class, same message, same guard. **[M]** After the fix, the deck writes it (see §7).

**Cantera YAML v1** does not export reactions at all in this path; **[R]**
`rmgpy/yaml_cantera1.py:137-145` only *refuses* charged species, and is disabled by default.

**RMS: no — it fails earlier, for an unrelated reason, and this ticket does not fix it.** **[R]**
`rmgpy/yaml_rms.py:144` builds the reactant list straight from `obj.reactants` with no electron
handling and no guard at all, so the electron question never arises. It never gets that far:
**[R]** `rmgpy/yaml_rms.py:146` serialises the kinetics first, and `:250-253` has no branch for
`VoronovEIArrhenius` or `BadnellRRArrhenius`. **[M]**

```
-- RMS (rmgpy/yaml_rms.py) --
   ValueError: Object of type <class 'rmgpy.kinetics.arrhenius.VoronovEIArrhenius'> does not
   have a defined conversion to ReactionMechanismSimulator format
```

**[M]** This is what killed the deck on the *second* save even after the Chemkin blocker was
fixed — the run got past model generation and died in `RMSWriter.update`. **This is reported as a
partial fix, not presented as done.** Two things follow, and both are stated rather than papered
over:

1. **The gap is a kinetics-coverage gap, not an electron-representation one.** It fires on the
   rate law before any equation is built, and it would fire identically on a reaction with no
   electron in it. Filling it means deciding how a Te-dependent rate law is represented on RMS's
   Julia side, which is a design question this ticket has no basis to answer — guessing one would
   be a correct implementation of a wrong design.
2. **The deck now turns the RMS writer off**, `generateRMSYAML=False` in
   `docs/i123-integration/input.py`, with the reason written next to it. That is a configuration
   statement of fact — RMS cannot represent this chemistry — not a way past a check. The writer is
   enabled by default **[R]** `rmgpy/rmg/input.py:1782` and would kill any plasma run using these
   rate laws. The same deck now turns the Cantera v2 writer **on**, so the second writer sharing
   `electron_balance` is exercised by the real mechanism rather than only by unit test.

---

## 6. The fix

`rmgpy/electron_balance.py`:

- **new** `get_placement_declaration(reaction)` — returns the owner's
  `(reactant_count, product_count)` or `None`. `None` means "no statement made"; a *malformed*
  declaration raises rather than falling back, because a silent fallback to the rule a declaration
  exists to override is how a wrong file gets written while the export reports success. Lazy
  import, so `electron_placement → electron_balance` stays a one-way dependency.
- **`expand_electrons`** — consults it. Declaration present: validate `reaction.electrons` against
  `product_count - reactant_count`, then place both sides from the declaration. Absent: the
  net-derived rule, byte-for-byte as before. `electrons == 0` still returns early, so a reaction
  already in explicit form — the reactor's view, or a mechanism read back from a file — passes
  through untouched and the helper is idempotent.
- module docstring records why the join is here.

`rmgpy/electron_placement.py`: `_place_declared_electrons`'s docstring corrected — it argued for
the split on a premise this ticket refuted, and leaving it would have the codebase contradicting
itself in the two files a reader would consult together.

Neither file is cythonized (**[M]** neither appears in `setup.py`'s `ext_modules`), so no `.so`
carries this logic; `make build` was run anyway and exits 0.

---

## 7. The written file, inspected

A green writer is not a correct file. **[M]** `<out>/chemkin/chem.inp` from the passing deck run:

```
ELEMENTS
	Ar  E  He  Li  N  Ne
END
SPECIES
    He  Ar  Ne  N2  e-(1)  Li(2)  [Lip](3)
END
...
[Lip](3)                E  -1Li  1          G    10.000  3000.000  794.00      1
...
REACTIONS    KCAL/MOLE   MOLES

Li(2)+e-(1)=>[Lip](3)+e-(1)+e-(1)                   1.490e+18 -0.267    162.150
    TDEP/e-(1)/   ! VoronovEIArrhenius exported as a modified-Arrhenius fit of k(Te) over
                    1.16e+04-2.321e+08 K; the original functional form is not representable
                    in Chemkin
END
```

Checked, item by item:

- **the electron is on both sides**, one incident and two liberated — the `(1, 2)` declaration,
  written out.
- **`E` balances**: left `Li(2)` E=0 + `e-(1)` E=+1 = **+1**; right `[Lip](3)` E=−1 + two
  electrons E=+1 each = **+1**. The species block agrees: `[Lip](3)` is declared `E -1 Li 1`.
- **the rate is not off by a factor of the electron density.** `A = 1.490e+18` is in
  `cm^3/(mol*s)` — second order — beside an equation with two reactants. Order of the coefficient
  and order of the equation agree, which is the entire content of the guard. Had it been written
  as `Li(2)=>[Lip](3)+e-(1)`, the same number beside one reactant would have overstated the rate
  by `1/[e] = 6.0 × 10⁶` at the deck's electron density (§2.3).
- **the loss is declared, not hidden.** **[M]** The exported Arrhenius reduction is
  `k(Te = 11 604.5 K) = 1.082e+08` against the Voronov coefficient's `1.293e+08`, a ratio of
  **0.837**. That 16% is the *reduction's* fit error over its own 11 604 K – 2.3 × 10⁸ K window —
  Chemkin has no Voronov form — and the file says so on the `TDEP` line beside the numbers. It is
  not a round-trip loss and not an order error.

**[M]** `<out>/cantera2/chem_annotated.yaml` from the same run, the native Cantera writer:

```yaml
reactions:
- equation: Li(2) + e-(1) => [Lip](3) + e-(1) + e-(1)
  type: two-temperature-plasma
  rate-constant: {A: 1490435700692.9358, b: -0.26677974840496976,
                  Ea-gas: 678436.9302780825, Ea-electron: 678436.9302780825}
```

Same placement, same order, and Cantera's native two-temperature form rather than a reduction.
`Ea-gas == Ea-electron` is deliberate and documented **[R]** `rmgpy/yaml_cantera2.py:801-806`: it
is the condition under which the two gas-temperature factors cancel and `k` is a pure function of
`Te`.

---

## 8. Round trip

**[M]** `test/rmgpy/i126ChemkinElectronOrderTest.py::TestChemkinRoundTrip::test_round_trip`
writes the mechanism, reads it back with `load_chemkin_file`, and shows one incident electron and
two product electrons survive, at two reactants, with `k` matching the written reduction to within
2% at 12 000 K, 20 000 K and 100 000 K (the Chemkin line carries four significant figures).

**The full trip is blocked one step short, by a defect this ticket did not introduce and does not
own.** RMG's Chemkin **reader** has no case for the `TDEP/<electron>/` auxiliary line RMG's own
Chemkin **writer** emits for every plasma rate law:

```
ChemkinError: 'e(1)' doesn't look like a collision efficiency for species TDEP in line
'TDEP/E(1)/'
```

**[M]** Established as pre-existing by reproducing it on **radiative recombination**, whose
declaration is one-sided and whose exported equation is therefore byte-identical before and after
this ticket's change — `Liplus(3)+e(1)=>Li(2)`. It is pinned as such by
`test_tdep_line_is_what_blocks_the_full_trip`, so it cannot regress quietly and cannot be
mistaken for this ticket's doing. The round-trip test drops that one auxiliary line — which is
what any Chemkin parser without plasma support sees — and round-trips everything else.

**[M]** The same gap appears in third-party Cantera: RMG's post-hoc `cantera_from_ck` step shells
out to `ck2yaml`, which fails on the same line with `could not convert string to float:
'e-(1)'`. Non-fatal — RMG logs it and continues, and the run still exits 0 — but it means the
`cantera_from_ck/` directory is empty for any plasma mechanism. The **native** Cantera v2 writer
(§7) is unaffected and is the one that produces a usable file.

---

## 9. Negative control

**[M]** `examples/rmg/minimal/input.py`, run twice on the same tree — once with
`rmgpy/electron_balance.py` at base `3294ed3fa` and once at the fixed version — into separate
output directories:

```
NEGCTL_BEFORE_EXIT=0
NEGCTL_AFTER_EXIT=0

IDENTICAL  chem.inp             (308924af2ca69954eddcba1853080886)
IDENTICAL  chem_annotated.inp   (d5108f38ff4634aec225776ae82c47ff)
IDENTICAL  chem_edge.inp        (6fc758854b2a3975c0c3d5216f7e8523)
IDENTICAL  species_dictionary.txt (c321b3b67de7ad4ac3bddfb23beea391)
```

Byte-identical, as the blast-radius argument in §4 predicts: an all-neutral mechanism never
reaches the declaration branch at all, because `electrons == 0` returns first.

**[M]** The unit suite, same scope the I-123 audit ran (`-m "not functional and not database"`):

```
= 2 failed, 2834 passed, 50 skipped, 131 deselected, 54 warnings in 158.48s =

FAILED test/rmgpy/electronPlacementTest.py::TestElectronPlacementResolver
       ::test_declaration_registry_is_explicit_and_closed
FAILED test/rmgpy/preflightDeckFamilyExclusionTest.py::PlasmaDeckFamilyExclusionTest
       ::test_plasma_deck_does_not_declare_family[docs/i102-quarantine/input.py]
```

**Both failures are the two the I-123 audit already recorded on this base**, by name and by
cause: §6.3 (`i119-rr-registry` shipped a fourth registry entry and updated two of the three
closed-registry assertions) and §6.2 (`i115-preflight-deck`'s deck sweep meets `i102-quarantine`'s
deliberate gate-demonstration deck). Neither is touched by this ticket, and no new failure
appears. The count rises from the audit's 2816 by the 17 tests added here.

---

## 10. Durable findings

1. **`Reaction.is_balanced` and the writers use two different rules for the electron**, and only
   the writers' rule can express ionisation. `is_balanced` counts the `e` pseudo-element by atom
   count **[R]** `rmgpy/reaction.py:1721-1723`; `get_species_electron_count` counts it as minus
   the net charge **[R]** `rmgpy/electron_balance.py:228-242`. This is *why* the canonical
   database form has no explicit electron — a fact that was recorded in the library's `longDesc`
   as a resolver constraint, when the binding constraint is actually the loader's balance check.
   Any future attempt at layer 2 hits this first.
2. **`KineticsLibrary.load` propagates `kinetics.electrons` onto the reaction unconditionally**
   **[R]** `rmgpy/data/kinetics/library.py:567-583`, so a library entry cannot opt out of the
   metadata representation even by writing `electrons = 0`.
3. **RMS export has no branch for either plasma rate law** **[R]** `rmgpy/yaml_rms.py:250-253`,
   and the RMS writer is on by default **[R]** `rmgpy/rmg/input.py:1782`. Any plasma deck that
   does not set `generateRMSYAML=False` dies on its first save, *after* model generation
   succeeds. This is the next blocker on this path.
4. **RMG's Chemkin reader cannot read RMG's Chemkin writer's `TDEP` line**, and neither can
   Cantera's `ck2yaml` (§8). Affects every plasma rate law, not just the ionisation shape.
5. **A plain `Reaction` has no `family` attribute at all** (`AttributeError` on assignment,
   since it is cythonized without one). Convenient here — it is what guarantees the reactor's
   placement view can never be re-expanded — but worth knowing before writing a test that sets
   `reaction.family`.

---

## 11. What I could not reach

1. **Whether the exported mechanism is *chemically* right past the equation.** The deck reaches
   `MODEL GENERATION COMPLETED` with 7 core species and 1 core reaction, and the recombination
   channel never enters the core (`This library reaction was not new`), so the written file
   contains the source and not the sink. Whether that is correct model behaviour or a
   flux/tolerance artefact is a model question, not an export question. **UNKNOWN.** Settled by:
   running the deck with the recombination reaction forced into the core and comparing the
   electron-density trajectory.
2. **Whether any *other* two-sided declaration will ever be added, and whether the validation in
   `expand_electrons` is right for it.** Today's registry has exactly one, so the two-sided path
   has one shape's worth of evidence behind it. The three-body recombination `(2, 1)` case is
   already named as blocked for an unrelated reason **[R]** `rmgpy/electron_placement.py:145-154` (unchanged by this ticket).
3. **Whether `Reaction.to_cantera` (the in-memory path, `reaction.py:312`) is now correct
   end-to-end.** It gets the right stoichiometry from the fixed converter, but **[M]** it then
   raises `NotImplementedError: Unable to set cantera kinetics for VoronovEIArrhenius` from the
   kinetics side, so the electron half cannot be observed through it. Not this ticket's gap; the
   same class of gap as RMS.
4. **The functional and database-marked suites.** Only `-m "not functional and not database"` was
   run, matching the audit's scope. **UNKNOWN**, and the audit's §10.4 already names this as a
   standing gap on this lineage.

---

## Reproducing this

```bash
conda activate rmg_env
cd /home/alon/Code/RMG-Py-i126-chemkin-electrons
export PYTHONPATH=/home/alon/Code/RMG-Py-i126-chemkin-electrons:$PYTHONPATH
python -c "import rmgpy; print(rmgpy.__file__)"                    # must be this worktree
python -c "from rmgpy import settings; print(settings['database.directory'])"
                                                                   # must be RMG-database-i123-integration/input
make build                                                         # never bare `make`

python -m pytest test/rmgpy/i126ChemkinElectronOrderTest.py        # 17 passed
python rmg.py -o <out> docs/i123-integration/input.py              # exit 0; §7 is <out>/chemkin/chem.inp
python rmg.py -o <out2> examples/rmg/minimal/input.py              # negative control, §9
python -m pytest -m "not functional and not database"              # §9: 2 failed, both pre-existing
```

Long-running commands were run as
`<command> > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)`; the RMS traceback in §5 is only
in `stderr.log`, not in `RMG.log`, which is why it is captured.
