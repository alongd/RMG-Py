# kMC rate-rule preparation from pinned training data

The compiler now prepares every loaded non-auto-generated family by calling
RMG's `add_rules_from_training(thermo_database=...)`, followed by
`fill_rules_by_averaging_up(verbose=True)`. Training preparation temporarily
clears species and polymer constraints and restores them even on failure.
An explicit `!training` request skips addition while retaining averaging.
Prepared families are reused once within a loaded database. Auto-generated
families are skipped entirely.

The artifact records the procedure, per-family rule counts, and database
commit `4a12d36fcdc193ede82c8d1ab5c1653495d445bc`. The changed compiler source
fingerprint changes artifact/cache identity. Affected family records store
RMG's template, full estimation comment, and the rule/training entries and
weights extracted by RMG itself. Detailed-balance reverse records retain
that forward source. No rate, thermo or library is chosen from a benchmark.

## Reproduction

All scratch, caches and logs are under `/home/alon/runs/i046-rules-from-training/`.
Use this checkout as the working directory, and the supplied environment:

```bash
export PYTHONPATH=$PWD
export PYTHONHASHSEED=0
export RMG_KMC_CACHE_ROOT=/home/alon/runs/i046-rules-from-training/cache
export RMG_KMC_ARTIFACT=$RMG_KMC_CACHE_ROOT/artifact/7491ed418f3ce2f0633688a9487cd12da278704b8f57e4f595176ddf80fb0296.json
export RMG_DATABASE_PATH=/home/alon/runs/i046-rules-from-training/database
export RMG_DATABASE_SHA=4a12d36fcdc193ede82c8d1ab5c1653495d445bc
export RMG_KMC_FAMILY_UNIVERSE=/home/alon/runs/i046-rules-from-training/family-universe.json
export MPLCONFIGDIR=/home/alon/runs/i046-rules-from-training/mpl
export PYTHONPYCACHEPREFIX=/home/alon/runs/i046-rules-from-training/pycache
export TMPDIR=/home/alon/runs/i046-rules-from-training/tmp
export COVERAGE_FILE=/home/alon/runs/i046-rules-from-training/.coverage-verified
```

Build: `/home/alon/anaconda3/envs/rmg_env/bin/python setup.py build_ext --inplace`.
Materialize only the existing I034 allowlist using `git show`, plus the family
name universe using `git ls-tree`:

```bash
/home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i046_probe/materialize_database.py \
  /home/alon/runs/i046-rules-from-training \
  > >(tee -a /home/alon/runs/i046-rules-from-training/materialize.stdout.log) \
  2> >(tee -a /home/alon/runs/i046-rules-from-training/materialize.stderr.log >&2)
```

The intended PS compilation command (do not rerun for verification):

```bash
PYTHONHASHSEED=0 /home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/compile_event_set_fixture.py \
  /home/alon/runs/i046-rules-from-training/database \
  /home/alon/runs/i046-rules-from-training/cache/artifact \
  --database-sha 4a12d36fcdc193ede82c8d1ab5c1653495d445bc \
  --family-universe /home/alon/runs/i046-rules-from-training/family-universe.json \
  > >(tee -a /home/alon/runs/i046-rules-from-training/compile.stdout.log) \
  2> >(tee -a /home/alon/runs/i046-rules-from-training/compile.stderr.log >&2)
```

`RMG_DATABASE_SHA` supplies the pin for atom-map checks on the materialized
snapshot, which has no Git metadata. The default suite
and selected direct-artifact slow tests use this file and launch no compile.
The new fast artifact regressions fail against the supplied old artifact.
The comparison probe reconstructs old rule estimates from the pinned family
files, checks their entire rate grids, and compares the new artifact by
chemical rewrite identity rather than provenance-dependent event IDs.
Nested inverse handles use the compiler's existing canonical link normalization;
artifact validation still checks the concrete partner links.
It reuses I044's benchmark, graph reconstruction and numerical checks.

```bash
/home/alon/anaconda3/envs/rmg_env/bin/python -m pytest test/rmgpy/kmc -q \
  -o cache_dir=/home/alon/runs/i046-rules-from-training/pytest-cache \
  --cov-report=html:/home/alon/runs/i046-rules-from-training/htmlcov-verified \
  > >(tee -a /home/alon/runs/i046-rules-from-training/suite-verified.stdout.log) \
  2> >(tee -a /home/alon/runs/i046-rules-from-training/suite-verified.stderr.log >&2)
RMG_KMC_SLOW=1 \
COVERAGE_FILE=/home/alon/runs/i046-rules-from-training/.coverage-slow-focused \
  /home/alon/anaconda3/envs/rmg_env/bin/python -m pytest \
  test/rmgpy/kmc/compilerRealTest.py -q \
  -k 'family_enumeration or j_ortho or exact_id_pack or c7_c9' \
  -o cache_dir=/home/alon/runs/i046-rules-from-training/pytest-cache \
  --cov-report=html:/home/alon/runs/i046-rules-from-training/slow-focused-htmlcov \
  > >(tee -a /home/alon/runs/i046-rules-from-training/slow-focused.stdout.log) \
  2> >(tee -a /home/alon/runs/i046-rules-from-training/slow-focused.stderr.log >&2)
/home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i046_probe/compare_artifacts.py "$RMG_KMC_ARTIFACT" \
  --database /home/alon/runs/i046-rules-from-training/database \
  --output /home/alon/runs/i046-rules-from-training/results.json --verify \
  > >(tee -a /home/alon/runs/i046-rules-from-training/verify-final.stdout.log) \
  2> >(tee -a /home/alon/runs/i046-rules-from-training/verify-final.stderr.log >&2)
/home/alon/anaconda3/envs/rmg_env/bin/python \
  test/rmgpy/kmc/fixtures/i046_probe/render_tables.py \
  /home/alon/runs/i046-rules-from-training/results.json \
  > >(tee -a /home/alon/runs/i046-rules-from-training/render-final.stdout.log) \
  2> >(tee -a /home/alon/runs/i046-rules-from-training/render-final.stderr.log >&2)
```

Initial result generation uses the same comparison command without `--verify`;
report generation uses the renderer's `--update`. Neither recompiles events.

## Per-family estimates and rate changes

Record counts include both directions. Default/root counts count directly
estimated forward records; detailed-balance inverses inherit their source
and are not a second independent estimate. “Default/root” means an original
rule whose description is `Default` (case insensitive), an estimate selecting
the full non-tree root template (including an averaged root), or an ATG
estimate from `Root`. Partial ancestor templates are not full roots.
The rate distribution includes both directions
of every matched record, with the same 600, 700 and 800 K grid nodes.
Some new estimates still use averaged ancestors. The representative migration
uses `RnH;C_rad_out_Cs2;XH_out` for template
`R4H;C_rad_out_Cs2;Cb_H_out`; its full RMG comment and contributors are
retained in `results.json` and the artifact.

<!-- BEGIN I046:families -->
| Family | Old records | New records | Old Default/root estimates | New Default/root estimates |
| --- | --- | --- | --- | --- |
| Disproportionation | 12680 | 12680 | 0 | 0 |
| H_Abstraction | 628 | 628 | 314 | 0 |
| R_Addition_MultipleBond | 598 | 598 | 299 | 0 |
| R_Recombination | 594 | 594 | 0 | 0 |
| intra_H_migration | 498 | 498 | 249 | 0 |
<!-- END I046:families -->

<!-- BEGIN I046:preparation -->
| Family | Rules loaded | Original Default entries | After training | After averaging |
| --- | --- | --- | --- | --- |
| H_Abstraction | 2 | 2 | 3121 | 9263 |
| R_Addition_MultipleBond | 1 | 1 | 2963 | 7392 |
| intra_H_migration | 12 | 1 | 451 | 6686 |
<!-- END I046:preparation -->

<!-- BEGIN I046:rate_changes -->
| Family | T (K) | Matched records | Median log10(new/old) | Minimum | Maximum |
| --- | --- | --- | --- | --- | --- |
| Disproportionation | 600.0 | 12680 | 0 | 0 | 0 |
| Disproportionation | 700.0 | 12680 | 0 | 0 | 0 |
| Disproportionation | 800.0 | 12680 | 0 | 0 | 0 |
| H_Abstraction | 600.0 | 628 | 1.90562 | -4.78415 | 9.49942 |
| H_Abstraction | 700.0 | 628 | 2.56444 | -2.95441 | 9.0823 |
| H_Abstraction | 800.0 | 628 | 3.08276 | -1.54834 | 8.80323 |
| R_Addition_MultipleBond | 600.0 | 598 | -6.03436 | -8.06314 | -0.396192 |
| R_Addition_MultipleBond | 700.0 | 598 | -5.40668 | -7.23095 | -0.25778 |
| R_Addition_MultipleBond | 800.0 | 598 | -4.9153 | -6.57974 | -0.141209 |
| R_Recombination | 600.0 | 594 | 0 | 0 | 0 |
| R_Recombination | 700.0 | 594 | 0 | 0 | 0 |
| R_Recombination | 800.0 | 594 | 0 | 0 | 0 |
| intra_H_migration | 600.0 | 498 | -5.18153 | -27.3182 | 7.48819 |
| intra_H_migration | 700.0 | 498 | -4.07872 | -23.4156 | 6.60839 |
| intra_H_migration | 800.0 | 498 | -3.26525 | -20.4886 | 5.96026 |
<!-- END I046:rate_changes -->

## Styrene propagation

The two channels are the same exact chemical rewrites selected by I044, now
matched by structure. Rates are the compiled tables in L/mol/s. IUPAC values
reuse the earlier PLP-SEC Arrhenius transcription; all temperatures here are
extrapolations beyond its measured range. The original primary radical proxy
and the bulk benzylic propagation benchmark remain different chemical/phase
targets. Ratios are a diagnostic comparison, not a tuning target or error bar.
RMG selects the exact `Cds-CbH_Cds-HH;CsJ-CsHH` rule, index 4324, rank 10,
generated from training reaction 1408 (`Aaron Vandeputte GAVs CBS-QB3`).
RMG's source extractor reports this as a rule contributor of weight 1;
its training origin is retained in the rule description and full comment.

<!-- BEGIN I046:propagation -->
| Proxy | T (K) | k_fwd (L/mol/s) | IUPAC extrapolation (L/mol/s) | Ratio |
| --- | --- | --- | --- | --- |
| end_radical | 600.0 | 2355.69 | 63257.2 | 0.0372398 |
| end_radical | 700.0 | 8875.13 | 160435 | 0.0553193 |
| end_radical | 800.0 | 25059.5 | 322433 | 0.07772 |
| end_radical@5 | 600.0 | 2355.69 | 63257.2 | 0.0372398 |
| end_radical@5 | 700.0 | 8875.13 | 160435 | 0.0553193 |
| end_radical@5 | 800.0 | 25059.5 | 322433 | 0.07772 |
<!-- END I046:propagation -->

<!-- BEGIN I046:sources -->
Proxy `end_radical`, template `Cds-CbH_Cds-HH;CsJ-CsHH`:

```text
From training reaction 1408 used for Cds-CbH_Cds-HH;CsJ-CsHH
Exact match found for rate rule [Cds-CbH_Cds-HH;CsJ-CsHH]
Euclidian distance = 0
family: R_Addition_MultipleBond
```

```json
{
  "entry": "Cds-CbH_Cds-HH;CsJ-CsHH",
  "exact": true,
  "rank": 10,
  "rules": [
    {
      "entry": {
        "index": 4324,
        "label": "Cds-CbH_Cds-HH;CsJ-CsHH",
        "rank": 10,
        "short_desc": "Rate rule generated from training reaction 1408. Aaron Vandeputte GAVs CBS-QB3"
      },
      "weight": 1
    }
  ],
  "training": []
}
```

Proxy `end_radical@5`, template `Cds-CbH_Cds-HH;CsJ-CsHH`:

```text
From training reaction 1408 used for Cds-CbH_Cds-HH;CsJ-CsHH
Exact match found for rate rule [Cds-CbH_Cds-HH;CsJ-CsHH]
Euclidian distance = 0
family: R_Addition_MultipleBond
```

```json
{
  "entry": "Cds-CbH_Cds-HH;CsJ-CsHH",
  "exact": true,
  "rank": 10,
  "rules": [
    {
      "entry": {
        "index": 4324,
        "label": "Cds-CbH_Cds-HH;CsJ-CsHH",
        "rank": 10,
        "short_desc": "Rate rule generated from training reaction 1408. Aaron Vandeputte GAVs CBS-QB3"
      },
      "weight": 1
    }
  ],
  "training": []
}
```
<!-- END I046:sources -->

## Ceiling, discovery and artifact identity

The probe independently re-estimates the selected graph thermochemistry,
checks every Kc grid node and detailed-balance reverse rate, and solves the
continuous 1 mol/L crossing. The continuous thermo ceiling and interpolated
artifact-grid crossing are different numerical quantities.

<!-- BEGIN I046:provenance -->
| Quantity | Reproduced value |
| --- | --- |
| New artifact SHA-256 | 7491ed418f3ce2f0633688a9487cd12da278704b8f57e4f595176ddf80fb0296 |
| New record count | 14998 |
| Tree records with identical rates/sources | 13274 |
| Old non-tree rate nodes reproduced | 24998 |
| Old artifact-grid ceiling (K) | 710.249 |
| New artifact-grid ceiling (K) | 710.249 |
| end_radical continuous thermo ceiling (K) | 710.020462 |
| end_radical@5 continuous thermo ceiling (K) | 710.020462 |
| Chemical records gained | 0 |
| Chemical records lost | 0 |
| Status changes | 0 |
| Old excluded channels | 2 |
| New excluded channels | 2 |
<!-- END I046:provenance -->

<!-- BEGIN I046:discovery -->
| Family | Old discovery records | New discovery records |
| --- | --- | --- |
| Disproportionation | 6340 | 6340 |
| H_Abstraction | 318 | 318 |
| R_Addition_MultipleBond | 300 | 300 |
| R_Recombination | 294 | 294 |
| intra_H_migration | 251 | 251 |
<!-- END I046:discovery -->

## Contract corrections and verification evidence

The contract's statement that all three non-tree rule files contain only
defaults is too broad: pinned `intra_H_migration/rules.py` contains one Default
rule (index 614) and 11 specific rules. For example, index 649 is labelled
`R2H_S;Cd_rad_out;Cs_H_out_OneDe` and described as
`Sumathy B3LYP/CCPVDZ calculations`; index 1018 is described as
`A. G. Vandeputte BMK/cbsb7 HO`. The probe saves the original entry metadata
in `results.json` and renders the original Default counts above. Actual old default
use is measured above rather than inferred from file size. Reaction discovery
uses graphs and family recipes; the compiler does not filter records by rate
thresholds.

Execution deviation: the default suite was started before the artifact was
available. An unmarked SSA test automatically launched another compiler
process on a cache miss. That process was stopped during junction-pair
reaction discovery, along with the premature suite (exit 143). It produced
no artifact. This breached the instruction to launch exactly one compile.
The intended compiler continued. All five artifact fixtures now accept
`RMG_KMC_ARTIFACT`; an invalid explicit file fails rather than launching a
fallback compile. The final suite is run only after the intended artifact
exists. The aborted attempt's logs remain under the task's cache at
`d027e21d0bb48bb0e624835815acfdd8ee96c3b4-4a12d36fcdc193ede82c8d1ab5c1653495d445bc-a6fe094da7fc97aad3ec3077d29a5f89c1cc6313f4f5446347722e97f2cd1fbb/`.

The intended compilation completed and produced the artifact above. No compile
was launched in the resumed thread. The initial resumed default-suite attempt
found that the atom-map SHA check expected Git metadata in the materialized
snapshot; `RMG_DATABASE_SHA` now supplies that pin, and the corrected run passed
the check. Its final logs are `suite-verified.stdout.log` and
`suite-verified.stderr.log`.

Both the initial comparison and its independent `--verify` rerun exited 0:

```text
I046 artifact hashes, compiler fingerprint and schema verified
I046 tree-family rate tables and sources: 13274 bit-identical records
I046 old non-tree rule rates reproduced; family distributions and discovery compared
I046 both propagation pairs: dimensional Kc, detailed balance and 710.02 K thermo ceiling verified
I046 three fresh public-RMG estimates: all 87 rate nodes and source comments verified
```

The report renderer exited 0 with `I046 all seven report blocks verified`.
The old-artifact negative control, selecting only the three representative
training-rule regressions, exited 1 as expected: `3 failed, 6 deselected in
13.90s`. Each failed on missing training-derived sources. This control loads
the old file directly and never compiles.

An optional exhaustive slow selection (`-k 'not frontier'`) was stopped with
exit 143 after its family-coverage and two ortho tests passed. Its archived
graph-invariance scan was still active after more than 15 minutes. The focused
selection above completes the direct-artifact checks that fit the dispatch's
optional scope. The comparison probe and fast artifact regression independently
cover all tree-family rate tables and sources. The broader all-record rewrite
and all-pair thermo slow checks were not completed. Independent generation
oracles and cross-process compile determinism were not rerun, and the frontier
test was excluded because it performs another compile.

The final default suite exited 0:

```text
177 passed, 48 skipped in 758.34s (0:12:38)
```

Its real evolving trajectory completed all 10,000 events across all five
families. The focused slow selection exited 0:

```text
5 passed, 7 deselected in 248.89s (0:04:08)
```

Those checks cover the complete pinned family universe, both ortho archive
checks, exact pack Kc and its corruption control, and C7/C9 inventory,
ceiling, inverse links and artifact corruption controls. All used the same
supplied artifact. No required work remains; no push was made.
