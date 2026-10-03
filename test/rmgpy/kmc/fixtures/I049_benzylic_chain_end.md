# Benzylic PS chain end: Phase 1 diagnosis and fix design

Branch: `i049-benzylic-chain-end`, starting commit `2e94caa6f`.
Database pin: `4a12d36fcdc193ede82c8d1ab5c1653495d445bc`.
Artifact: `/home/alon/runs/i046-rules-from-training/cache/artifact/7491ed418f3ce2f0633688a9487cd12da278704b8f57e4f595176ddf80fb0296.json`.
Scratch: `/home/alon/runs/i049-benzylic-chain-end/`.

This is a discovery defect, accompanied by an independent catalogue defect.
The normal head-to-tail propagation and its exact unzipping inverse are absent.
Other terminal benzylic radicals and their inverse reactions do exist. No
production code, event set, rate selection, or thermo selection is changed here.

## Reproduction / Verifier

From `/home/alon/Code/RMG-Py-kmc-i049-benzylic-chain-end`:

```bash
export PYTHONPATH=$PWD
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONHASHSEED=0
export RMG_DATABASE_PATH=/home/alon/runs/i046-rules-from-training/database
export RMG_DATABASE_SHA=4a12d36fcdc193ede82c8d1ab5c1653495d445bc
export RMG_KMC_ARTIFACT=/home/alon/runs/i046-rules-from-training/cache/artifact/7491ed418f3ce2f0633688a9487cd12da278704b8f57e4f595176ddf80fb0296.json
export RMG_KMC_SLOW=0
export MPLCONFIGDIR=/tmp/i049-mpl
export COVERAGE_FILE=/tmp/i049.coverage
P=/home/alon/anaconda3/envs/rmg_env/bin/python
S=/home/alon/runs/i049-benzylic-chain-end
ulimit -v 8388608

cp /home/alon/Code/RMG-Py-kmc-i046-rules-from-training/rmgpy/solver/settings.pxi rmgpy/solver/settings.pxi
$P setup.py build_ext --inplace > >(tee -a "$S/build.stdout.log") 2> >(tee -a "$S/build.stderr.log" >&2)
$P test/rmgpy/kmc/fixtures/i049_probe/direct.py > >(tee -a "$S/direct.stdout.log") 2> >(tee -a "$S/direct.stderr.log" >&2)
$P test/rmgpy/kmc/fixtures/i049_probe/artifact.py --output "$S" > >(tee -a "$S/artifact.stdout.log") 2> >(tee -a "$S/artifact.stderr.log" >&2)
$P test/rmgpy/kmc/fixtures/i049_probe/runtime.py "$S" > >(tee -a "$S/runtime.stdout.log") 2> >(tee -a "$S/runtime.stderr.log" >&2)
$P test/rmgpy/kmc/fixtures/i049_probe/rates.py --output "$S/rates.json" > >(tee -a "$S/rates.stdout.log") 2> >(tee -a "$S/rates.stderr.log" >&2)
$P -m pytest test/rmgpy/kmc/compilerTest.py test/rmgpy/kmc/rulesFromTrainingTest.py \
  test/rmgpy/kmc/stateTest.py -k 'not compiled_representative and not tree_family and (not stateTest or terminal_backbone_form or evolving_coproduct_radical or evolving_recount)' \
  -q -o addopts='' -o cache_dir=/tmp/i049-pytest-cache \
  > >(tee -a "$S/tests.stdout.log") 2> >(tee -a "$S/tests.stderr.log" >&2)
```

The supplied starting probe was also rerun, in its required directory:

```bash
cd /home/alon/runs/i046-benzylic-probe
/home/alon/anaconda3/envs/rmg_env/bin/python corr.py \
  /home/alon/runs/i046-rules-from-training/cache/artifact/7491ed418f3ce2f0633688a9487cd12da278704b8f57e4f595176ddf80fb0296.json \
  > >(tee -a /home/alon/runs/i049-benzylic-chain-end/corr.stdout.log) \
  2> >(tee -a /home/alon/runs/i049-benzylic-chain-end/corr.stderr.log >&2)
```

Each artifact process loads that JSON once. `runtime.py` reads only the small
`runtime-records.json` subset written by `artifact.py`. `rates.py` checks every
materialized Python/text file in the family/thermo paths it loads against `git show`
at the pin, loads only `R_Addition_MultipleBond`, and makes two product-filtered
family-generation calls. None of these commands compile an event set. The
pytest selection explicitly excludes the cached old-artifact comparisons and
all expensive real-corpus state tests.

## Direct species and site-type diagnosis

`rmgpy/kmc/compiler.py:1133` puts `radical_end=True` only on the first carbon,
`[CH2]`; the last carbon is always the closed-shell benzylic `C(Ph)` unless
`radical_unit` makes it `[C](Ph)`. That latter option describes an internal
tertiary radical but makes a hydrogen-free terminal radical at the last unit,
so it is not a substitute for terminal `CH•(Ph)`.

`ps_proxy_set` builds `end_radical` and `end_radical_short` from that primary
head (`compiler.py:718`). Its declarations include only that end in the
unimolecular, end/end, junction/end, end/styrene and end/pristine inputs
(`compiler.py:755`, `compiler.py:769`, `compiler.py:776`, `compiler.py:783`,
`compiler.py:793`). The base size-3 styrene input uses a size-2 radical; the
extended input uses size 5 and therefore has a size-6 addition product
(`compiler.py:784`, `compiler.py:836`). All 22 direct proxy inputs lack a
terminal benzylic radical. This confirms cause 1.

The omitted ordinary single-radical species include:

| Repeat units | Missing benzylic-tail input SMILES | Formula |
| --- | --- | --- |
| 1 | `C[CH](c1ccccc1)` | C8H9 |
| 2 | `CC(c1ccccc1)C[CH](c1ccccc1)` | C16H17 |
| 3 | `CC(c1ccccc1)CC(c1ccccc1)C[CH](c1ccccc1)` | C24H25 |

The size-2 input plus styrene should give the size-3 member of this same
series; the size-3 input gives size 4. The missing discovery declarations are
the benzylic analogues of the above end types, including its shorter builder
and mixed benzylic/primary termination pair. Exact names proposed below.

The catalogue independently identifies `(site tuple)` with its index reversal
at `compiler.py:1164`. On this chain that reversal exchanges a carbon bonded
to phenyl with one having no phenyl neighbor. It is not a graph automorphism.
`direct.py` enumerates every 0-, 1- and 2-radical hydrogen-deletion candidate
and compares actual RMG graphs by `is_isomorphic`, without index reflection.
All 37 candidates are mutually nonisomorphic. The catalogue retains only 23:

| n | Stored | Correct | Missing single-radical site indices | Missing double-radical site tuples |
| --- | --- | --- | --- | --- |
| 1 | 3 | 4 | 1 (benzylic tail) | none |
| 2 | 7 | 11 | 2 (interior CH2), 3 (benzylic tail) | (1,3), (2,3) |
| 3 | 13 | 22 | 3 (interior benzylic), 4 (interior CH2), 5 (benzylic tail) | (1,5), (2,4), (2,5), (3,4), (3,5), (4,5) |

The direct log prints the exact SMILES, radical-site tuples and formulas for
every stored and missing graph, rather than just the counts. This confirms
cause 2. Its causal scope is narrower than the dispatch suggested: production
calls the catalogue only when serializing artifact metadata at
`compiler.py:2383`; discovery uses `self.proxies` and family generation
(`compiler.py:2188`). Catalogue validation only checks size and the termination
bound (`compiler.py:1369`). Fixing catalogue dedup alone adds **zero reaction
records** and does not repair discovery.

## What the artifact actually contains

The artifact has 14,998 enabled records, 7,503 discovery entries, and two
primary-end ceiling pairs. The supplied `corr.py` output was reproduced.
The independent structural probe distinguishes terminal benzylic radicals
(one H, one phenyl carbon neighbor, one non-ring carbon neighbor) from the
abundant internal tertiary benzylic radicals.

| Family | All records | Records containing a terminal benzylic product | Records containing a terminal benzylic reactant |
| --- | --- | --- | --- |
| R_Addition_MultipleBond | 598 | 4 | 4 |
| Disproportionation | 12,680 | 0 | 0 |
| H_Abstraction | 628 | 0 | 0 |
| R_Recombination | 594 | 22 | 22 |
| intra_H_migration | 498 | 3 | 3 |

These columns count records containing the structure, not net radical creation;
linked directions can appear in both columns. The standard tail series above
has zero occurrences on either side at n=1,2,3 and 5. There are **zero**
terminal-benzylic-radical + styrene addition records in the whole artifact.
The inverses of normal head-to-tail additions are consequently absent too.

Two existing primary-end additions make a terminal benzylic product, and
their two inverses release styrene. The short one is:

```text
[CH2]C(CCc1ccccc1)c1ccccc1 + C=Cc1ccccc1
  <=> [CH](CCC(CCc1ccccc1)c1ccccc1)c1ccccc1
template: Cds-HH_Cds-CbH;CsJ-CsHH
forward: evt_fa525f376be00ccd5d79c314618e6a408c783e90cdc07e1f8fd6a9f2f23e8629
inverse: evt_675c666a0e6b96da425c7469c054c918b8883c19b11ce61bc2f20bfc8aefe4ad
```

Here styrene adds to the primary radical's head, placing two unsubstituted
backbone carbons consecutively. The benzylic product is real but does not have
the alternating head-to-tail backbone. It cannot validate ordinary PS
propagation. The existing ceiling pair instead selects addition at styrene's
substituted carbon, retaining a primary radical. `artifact-results.json`
contains both exact ceiling rewrites and their IDs and temperatures.

## Runtime representation, rewrites and grammar

`Strand` stores position-dependent features, atom references, and an atom
graph (`state.py:67`); the parser retains radical, charge and bond-order data
(`state.py:170`). There is no restriction to a primary radical or to one
backbone carbon parity. A benzylic radical at either terminal backbone position
is representable without a schema change.

`SiteIndex` matches record graphs against current components (`ssa.py:596`)
and obtains site-type strings from each record (`ssa.py:434`, `ssa.py:648`).
It does not chemically convert an end based on its name. Preflight checks
site-type strings and matches the full chemical graph (`state.py:550`,
`state.py:630`). `radical_site_class` calls a radical at position 0 or length−1
an `end`, independently of primary/benzylic chemistry (`ssa.py:1149`). The
MET population keeps `(site_type, end_or_mid)` bins (`ssa.py:1423`); graph
rewrites set the mapped radical atom directly (`met.py:883`).

Compiler reverse-site labeling currently calls ring radicals `junction_radical`,
non-ring radicals with at least three C neighbors `interior_radical`, and
every remaining monoradical `end_radical` (`compiler.py:1753`). Thus a
terminal benzylic radical receives the generic `end_radical` label, and an
internal CH2 radical may receive that same label. This conflation affects
labels, not chemical graph storage; runtime's positional end/mid test remains
separate. Introducing a first-class benzylic label requires updating this
classifier as well as proxy declarations.

Actual atom edits retain their mapped position/chemistry (`state.py:858`).
Backbone bond cleavage derives a strand cut (`state.py:832`); splitting copies
features/atom graphs and repositions UUID references (`state.py:983`). Joining
terminal strands can reverse positions to orient the join, but preserves the
atom graph and radical features (`state.py:900`, `state.py:932`). It does not
exchange phenyl/nonphenyl atoms or turn a benzylic radical into a primary one.

The compiled routes that can currently supply terminal benzylic ends are:

* **Scission:** inverse `R_Addition_MultipleBond` can release a non-styrene
  alkene and a benzylic-ended fragment; inverse `R_Recombination` can cleave
  a backbone bond into radical fragments. These are not the missing normal
  styrene-unzipping channel. Examples and actual execution are recorded below.
* **Recombination:** forward radical recombination consumes the reacting
  radicals; its inverse homolysis supplies radicals, including a benzylic
  end where the mapped cleavage exposes one. Other spectator radicals, if
  present, are retained. It does not generate a surviving radical by ordinary
  two-monoradical combination.
* **Disproportionation:** forward termination consumes the two reacting
  radicals, so it does not create a surviving benzylic end. Its inverse could
  in principle create such an end from the right saturated/unsaturated graphs,
  but the compiled inventory has zero terminal-benzylic participants in this
  family. It must not be claimed as an existing source here.
* **H migration:** three product-containing records have terminal benzylic
  radicals, with nonstandard ring-localized/featured contexts. H migration
  moves the radical via mapped atom edits; no ordinary tail seed is present.
* **H abstraction:** in principle removing one H from the closed-shell
  benzylic tail would give the desired seed, but this artifact contains zero
  such participants. Adding a benzylic proxy must allow its inverse/transfer
  chemistry to be discovered rather than assuming the existing declarations
  already supply it.
* **Primary-end addition:** the channel explicitly shown above supplies a
  benzylic radical with a backbone defect.

For linked reversible records, `state.py:1269` bypasses coproduct sink debit;
the radical fragments remain active state graph nodes. Unlinked coproducts
can be intentionally debited into `sink_radicals` and their graph radical
removed (`state.py:1295`). That is explicit coproduct accounting, not a
primary/benzylic conversion. The runtime probe checks graph/ledger radicals,
inverse eligibility and chemical round trips on the linked examples.

`runtime.py` reproduced these four actual executions, each followed by its
compiled inverse. Graph and ledger radical totals agreed throughout;
`sink_radicals` remained zero, and every inverse restored the initial chemical
SMILES multiset and radical count. Cut offsets are zero-based backbone positions.

| Producer event ID | Channel | Strand lengths before → after | Radicals before → after | Derived cut |
| --- | --- | --- | --- | --- |
| `evt_0396caccabe4f92cb2c8878c5bcddbfe99d9c9f8068fdfa7ab11fd7738cdb4ae` | recombination inverse / homolysis | [6] → [3,3] | 1 → 3 | 2 |
| `evt_fa525f376be00ccd5d79c314618e6a408c783e90cdc07e1f8fd6a9f2f23e8629` | primary-end addition | [2,4] → [6] | 1 → 1 | none |
| `evt_1accc40f8f0d48777caf26d54d237f17b302676e5952bb781473580bc9152ac0` | featured H migration | [6] → [6] | 2 → 2 | none |
| `evt_d0ccbc31534990552476f319684f33f65b04e7d6b8896723621b6b749403faee` | featured addition inverse / scission | [10] → [5,5] | 2 → 2 | 4 |

All four after-states retain an actual terminal benzylic radical. These are
small exact-record checks, not a claim that every compiled reaction or a long
SSA trajectory was tested. The independent n=1–3 tail seeds each yielded one
site match with the new proposed label and an `end` positional class.

There is no implemented whole-component G8 chain-grammar recognizer in
`rmgpy/kmc/`. The older fixture identifies this as a gap:
`test/rmgpy/kmc/fixtures/i028_registered_verifier_mapping.md:13` and `:23`.
The current code search still finds no `accepts_g8`, `grammar_class` or G8
terminal schema in production. Therefore **atom-state representability is
supported; a formal grammar-closure guarantee is not implemented and cannot
be certified by this Phase 1 probe**. Graph matching/rewrite preservation
should be the immediate runtime regression, with any later formal grammar
work scoped separately.

## Design-intent search

Searched production kMC modules, the independent completeness oracle, all
fixture Markdown, and repository design Markdown under `docs/`, plus the
documentation/README text, for end/benzylic/head-to-tail/mirror/omission terms.
No rationale for excluding the benzylic end was found in that scope. This is
a statement about those inspected files, not inaccessible campaign history.

The relevant prior acknowledgement is
`test/rmgpy/kmc/fixtures/I039_tc_gap_probe.md:287`: its terminal-benzylic control
“is not compiled into the event set”. It is described as a diagnostic control,
not a chemical reason to exclude it. `I044_kp_benchmark.md:118` explicitly
distinguishes the terminal-primary compiled proxy from normal benzylic
propagation. Repository design guidance points the other way:
`docs/proxy_reaction_reality_rules.md:43` says “Caps are *allowed*”, and
`docs/multi_pool_design.md:322` includes chain-end unzipping.

The independent completeness oracle repeats the primary-only end input
(`test/rmgpy/kmc/completeness_oracle.py:32`) and its original scheduling
(`:43`). This explains how a chemically incomplete proxy inventory can pass
an independently implemented generation comparison: both enumerate the same
omitted class. The catalogue size test checks accounting, not coverage of
all distinct structures (`test/rmgpy/kmc/compilerTest.py:371`).

Search commands (no database/data tree searches):

```bash
rg -n 'linear_ps|short_ps|end_radical|ceiling_pairs' rmgpy/kmc/compiler.py
rg -n 'benzylic|head.to.tail|omit|exclude|palindrom|mirror|chain.end' \
  rmgpy/kmc test/rmgpy/kmc/fixtures --glob '*.md' --glob '*.py'
rg -n 'benzylic|head.to.tail|end_radical|chain.end|kMC|kmc' \
  docs documentation/source README* --glob '*.md' --glob '*.rst'
rg -n 'G8|accepts_g8|GRAMMAR_TERMINAL_SCHEMA|grammar_class' rmgpy/kmc
rg -n 'short_ps_molecule_catalogue|short_molecule_catalogue' rmgpy/kmc
```

The G8 search returns no production matches (rg exit 1).

## Minimal correct Phase 2 change set

1. Add a separate `benzylic_ps_end_smiles(n)` fragment helper, preserving
   existing `radical_end=True` primary-head behavior. The new helper
   must emit `CC(Ph)` repeated n−1 times plus `C[CH](Ph)`. Do not implement
   it by blindly reflecting carbon indices or by terminal `radical_unit`.
   Validate incompatible multiple features and end choices.
2. Retain all existing primary-end proxies. Add `benzylic_end_radical`,
   `benzylic_end_radical_short`, and the six discovery declarations
   `benzylic_end_radical`, `benzylic_end_radical+styrene`,
   `benzylic_end_radical+pristine`, `benzylic_end_radical+end_radical`,
   `benzylic_end_radical+benzylic_end_radical`, and
   `junction_radical+benzylic_end_radical`, with size-3 and size-5 contexts.
   Keep all five candidate families on the new inputs; do not rank chemistry
   using PLP-SEC agreement. Update `_direction_proxy` to distinguish terminal
   benzylic from primary and interior non-ring radicals from connectivity.
   This is the end-focused minimum; existing omissions such as mixed
   interior/end termination remain a separate completeness question.
3. For the base styrene declaration retain the short-to-long n=2→3 convention.
   Add a direct n=3→4 regression and a matching shorter extended benzylic
   input n=4→5 so **both** benzylic ceiling anchors actually have declared
   product proxies. The current primary extended styrene input n=5→6 may
   remain its own channel. Test exact reactant/product graphs, not a site-type
   substring, and explicitly locate the reacting benzylic radical and styrene
   carbon. Ordinary head-to-tail growth keeps a terminal `CH•(Ph)` radical
   and the alternating backbone in the product.
4. Remove index-reflection dedup from the catalogue. Enumerate the small
   bounded combinations and deduplicate actual featured molecule graphs by
   canonical structural key plus RMG isomorphism confirmation, or simply
   retain every combination for this PS series, for which the direct probe
   proves all are distinct at n≤3. A future generalization may quotient only
   by proven automorphisms of the substituted base molecule. Keep feature
   counts, formula, deterministic ordering and the same 37 termination bound;
   radius-1 size becomes 37. The catalogue correction alone costs 14 metadata
   entries, no event records.
5. Make `ps_ceiling_pairs` select **head-to-tail benzylic** propagation and
   its exact linked inverse. Preserve the two primary-end pairs under an
   separate `ps_primary_end_ceiling_pairs` field. Set the legacy scalar
   `ps_ceiling_temperature_K` from an explicitly designated benzylic anchor,
   rather than the minimum over unrelated channels. Update downstream
   fixture/probe consumers that currently assume two primary entries. This
   additive field is the proposed Phase 2 schema change.
   No new rate model or thermo selection is required. Continue estimating one
   family direction and deriving the exact inverse using the existing
   reference-state `Kc` (`compiler.py:1796`, `compiler.py:1901`). Verify
   `kf/kr=Kc` over the grid and dimensional units for each chemical pair.
6. Add independent fast coverage tests for both end orientations, n=1–3
   exhaustive catalogue graph coverage, n=2→3 and n=3→4 head-to-tail maps,
   regioisomer separation, reverse-site labels, and seed/site-index/end-class
   behavior. Extend `completeness_oracle.py` with independently written
   benzylic fragments and the new pair declarations, then rerun its real
   comparison only when the manager schedules the Phase 2 compile. Runtime
   regressions should apply a known mapped head-to-tail addition/unzip pair
   and prove UUID, ledger, radical and alternating-backbone preservation.

Adding only `+styrene` and omitting benzylic uni/transfer/termination inputs
would restore one growth channel but strand the chain-carrying radical's
other reactions. Conversely catalogue-only correction cannot restore growth.

## Recompile cost estimate, not a reproduced compile

Current records by family are in the artifact table above; disproportionation
accounts for 84.55% of the inventory. The comparable existing primary-end
declarations, summed over both contexts, give this scheduling estimate:

| New declaration analogue | Existing records used as weight |
| --- | --- |
| benzylic unimolecular ← end_radical | 180 |
| benzylic + styrene ← end_radical + styrene | 742 |
| benzylic + pristine ← end_radical + pristine | 2,860 |
| benzylic + primary ← end_radical + end_radical | 1,508 |
| benzylic + benzylic ← end_radical + end_radical | 1,508 |
| junction + benzylic ← junction_radical + end_radical | 2,859 |
| Total additional record estimate | 9,657 |

This yields about 24,655 total records and 34 declarations versus 22 now.
It is an analogy, not measured generation: resonance localization, degeneracy,
regioisomers, global dedup and the shorter extended styrene input can change
the result. A deliberately broad 0.5–1.5 multiplier on the added inventory
gives about **4,800–14,500 additional records**, or **19,800–29,500 total**.
The two core head-to-tail additions used for ceiling anchors alone require
at least two forward/inverse pairs (four records) before duplicate collapse.

The serialized artifact is 264,805,570 bytes (252.54 MiB), averaging about
17,656 bytes per record including shared metadata. Linear scaling gives
roughly 435 MB (415 MiB) centrally and 350–521 MB across that range. These
are disk-size estimates. The artifact structural probe measured about
934 MiB peak RSS while holding parsed records and distinct graph caches;
linear reader scaling suggests about 1.5 GiB centrally. **Compiler peak RAM
cannot be obtained from artifact composition**: it also retains RMG reaction
objects, mapping/thermo caches and artifact copies. The original artifact plus
two deep copies at `compiler.py:2455` would alone budget about 4.6 GiB at the
central reader scaling, before those RMG objects/caches. Reserve roughly 4–8 GiB
for a serial compile as a provisional scheduling envelope, measure RSS and
stop/replan if it reaches the cap. This is not an assurance that 8 GiB suffices.

Likewise neither the artifact nor its compile stdout records elapsed time or
peak RSS. The dispatch's eight-hour compile is a planning baseline, not a
timing measurement of this artifact. If that baseline is comparable, scaling
record work by 1.64 predicts about **13 h**; 1.32–1.97 gives **10.6–15.7 h**.
Rate-rule preparation is largely fixed, while molecular mapping and resonance
can dominate nonlinearly; allow **10–24 h**, at most four cores and 8 GiB, when
the manager schedules a single compile. No real event-set compile was performed
in Phase 1; the selected fast compiler tests use small mock inventories.

## Prepared-database head-to-tail rate pre-estimate

Both fragment calls found one merged reaction to the explicitly specified
head-to-tail product. Both selected template
`Cds-HH_Cds-CbH;CsJ-CbCsH`, a forward **rate-rule estimate**, with degeneracy
1.0, exact rule match and Euclidean distance 0. The restored rule is index
3253, rank 10, weight 1. Its comment identifies training reaction **337**:

```text
From training reaction 337 used for Cds-HH_Cds-CbH;CsJ-CbCsH
Exact match found for rate rule [Cds-HH_Cds-CbH;CsJ-CbCsH]
Euclidian distance = 0
family: R_Addition_MultipleBond
```

The original pinned entry is
`input/kinetics/families/R_Addition_MultipleBond/training/reactions.py:6240`
in the materialized database and at the same `git show` pin. It is labeled
`C8H8 + C8H9 <=> C16H17`, rank 10, description “Aaron Vandeputte GAVs CBS-QB3”.
Its long description at `:6255` says it was **converted to a training
reaction from this rate rule**. Thus it is a restored estimated rule,
not a direct PLP-SEC experiment or independent long-chain validation.
The script captures the original depository entry as `training_origin`;
RMG's source extractor itself lists a rule and an empty `training` array
for this comment format. The origin is not silently inferred from that array.

In both cases RMG returns
`ArrheniusEP(A=(926,'cm^3/(mol*s)'), n=2.41, alpha=0, E0=(31547.4,'J/mol'))`,
with Tmin=300 K and Tmax=1500 K. In L/mol/s the same expression is

```text
k(T) = 0.926 * (T / 1 K)^2.41 * exp(-31547.4 / (R*T))  L/mol/s
```

The call to `kinetics.get_rate_coefficient(T)` exactly matches the compiler
at `rmgpy/kmc/compiler.py:1556`, followed by the SI→L factor 1000. No
enthalpy conversion, benchmark fitting or alternative source is inserted;
alpha=0 makes the returned estimate independent of reaction enthalpy.
Preparation reproduced 1 loaded rule → 2,963 after training → 7,392 after
averaging, using `prepare_rate_rules` and the same primary thermo library as
the compiler (`compiler.py:136`, `test/rmgpy/kmc/compile_event_set_fixture.py:79`).
The final rate command verified 20 pinned files and completed in 25.31 s.

| Fragment reactant → product units | T (K) | RMG k (L/mol/s) | Given IUPAC Arrhenius k (L/mol/s) | RMG / IUPAC |
| --- | --- | --- | --- | --- |
| 2 → 3 | 600 | 8,233.170073 | 63,257.201921 | 0.130153877 |
| 2 → 3 | 700 | 29,461.174176 | 160,434.678657 | 0.183633454 |
| 2 → 3 | 800 | 80,032.529857 | 322,432.976207 | 0.248214469 |
| 3 → 4 | 600 | 8,233.170073 | 63,257.201921 | 0.130153877 |
| 3 → 4 | 700 | 29,461.174176 | 160,434.678657 | 0.183633454 |
| 3 → 4 | 800 | 80,032.529857 | 322,432.976207 | 0.248214469 |

The IUPAC column uses exactly the dispatch's A=4.27e7 L/mol/s and
Ea=32.5 kJ/mol with `rmgpy.constants.R`. It is only a numerical diagnostic.
The estimates happen to be identical because both fragments select the same
local rule and degeneracy. That does not establish real chain-length
independence or equality of a gas-fragment estimate and bulk PLP-SEC kinetics.
No thermo ceiling or reverse rate for these **new** reactions was computed
in this phase; detailed balance remains the stated Phase 2 verifier.

## Evidence and limits

The requested native build exited 0. The supplied `corr.py` exited 0 and
reproduced 14,998 enabled records, two primary-only ceiling pairs, and only
the two primary-to-benzylic styrene-addition channels plus their inverses.
The committed direct, artifact and runtime commands exited 0 with:

```text
PASS: n=1–3 catalogue 23/37, 14 distinct graphs omitted; 22 proxies omit benzylic tail
PASS: zero benzylic-end+styrene additions; existing ceiling pairs use primary ends
PASS: three benzylic-end seeds represented/indexed; 4 real producer/inverse round trips
PASS: two head-to-tail reactions estimated after pinned training preparation
23 passed, 24 deselected in 4.21s
```

The pytest result is from exactly the selection listed above; no event-set
compile and no `RMG_KMC_SLOW=1` suite was launched. The structural artifact
probe took 9.23 s with 956,996 KiB peak RSS in its final run. Full results and
logs are under the scratch path. The native build initially lacked importable
extensions; direct verification passed after the local products were built.
The first rate-probe product assertion compared aromatic family graphs with
raw Kekulé builder graphs and was corrected to compare canonical aromatic
SMILES. Its successful output is recorded in the rate section.

The contract's two suspected problems are both defects, but they are not
two independent eliminations in the discovery pipeline: catalogue dedup only
affects serialized coverage metadata. “No benzylic chain end” is too broad
for the artifact/runtime; the absent class is normal head-to-tail growth
from the benzylic end, and the standard H-capped tail seeds. Existing
benzylic-ended scission products, defective-backbone addition products and
featured migration products are present. The runtime's formal “grammar”
is an unimplemented recognizer, although the actual atom state can represent
the end. Compile wall time and peak RAM are projections, not reproduced
measurements. Phase 2 implementation and recompilation remain manager-gated.
