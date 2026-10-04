# PLP-SEC head-to-tail propagation library

## Implementation and scope

Implementation commit: `0da379aa5497ca6a03c86abc7ffafbd5009a9a14`.

`rmgpy/kmc/kinetics_library.py` carries one Python data mapping and a loader
that returns a deep copy. The entry contains the citation, units, temperature
range, extrapolation note, matching scope, and sensitivity alternative. The
compiler fingerprints the module alongside its existing source files.

The owner-selected expression is A=4.27e7 L/(mol*s), n=0, Ea=32.5 kJ/mol.
The original [Buback et al., Macromol. Chem. Phys. 196, 3267–3280 (1995)](https://doi.org/10.1002/macp.1995.021961016) reports −12 to 93 °C
(261.15–366.15 K), Eq. (5b), p. 3276. The
[publisher PDF in the authors' institutional repository](https://pure.tue.nl/ws/portalfiles/portal/1333387/617494.pdf) was checked.
The entry uses this original window, rather than the −12 to 90 °C window in
later summaries. Use at 600–800 K is explicitly an Arrhenius extrapolation.
The paper measures a bulk, low-conversion benchmark and advises avoiding
very short chains in the PLP analysis (p. 3269); it does not separately measure
the two- and four-unit proxy rate coefficients. Application to these finite
proxies is the owner's specified rate-model choice.

The matcher checks full chemical graphs **and** mapped edits, never templates
or declared site strings alone. It requires unsubstituted styrene plus an
ordinary H-capped alternating PS chain with at least two repeat units,
CH3–CH(Ph)–…–CH2–CH•(Ph). The product must be the same intact chain extended
by CH2–CH•(Ph). Isomorphism uses aromatic reference graphs. Exactly five
operations must occur: replace styrene's CH2=CHPh double bond with a single
bond, form the old radical–styrene CH2 single bond, remove the old radical,
and place the new radical on styrene CHPh. No other edit is allowed.

Independent mapped fixtures include n=2,3,4,6 growth with reversed reactant
ordering and misleading template/site labels, primary-end and head-to-head
look-alikes, actual ring addition, and all four other candidate families.
For excluded reactions both rate tables and rate sources are compared with
the disabled compiler. A correct product with a corrupted radical edit is
also rejected. The production declaration inventory's ordinary n=2→3 and
n=4→5 channels are the two expected library matches after recompilation.

## Precedence, inverse, and sensitivity

`EventSetCompiler(..., use_plpsec_library=False)` disables the library.
With no explicit constructor choice, `RMG_KMC_PLPSEC_LIBRARY=0` disables it;
`1` enables it and is the default. Invalid environment values fail explicitly.
The constructor takes precedence over the environment. The existing fixture
compiler therefore needs no modification to select a later sensitivity arm.

The RMG estimate is still evaluated and its template, source extraction,
original direction/units, and complete grid table are retained under
`rate_source.replaced_rmg_estimate`. `propagation_k_table` records that
alternative normalized to addition if RMG estimated unzipping. The library
is installed once per matched physical propagation reaction, with A converted
to 4.27e4 m^3/(mol*s); the measured coefficient receives no additional path
degeneracy factor. Original degeneracy metadata remains intact.

Kc comes from the existing reference-thermo provider. The inverse remains
kf/Kc, in s^-1, and includes the forward library provenance. If the family
estimated the reverse direction, the compiler reverses that orientation and
Kc before installing propagation. Tests cover both estimation directions.
Both on and off arms have identical mock ceiling temperatures, including the
designated anchor. No real new ceiling temperature is claimed.

Every matched forward source carries the library entry, citation, and replaced
estimate; its inverse carries that forward source. All records and the artifact
share `provenance.kinetics_libraries.styrene_plpsec`, including enabled state,
entry, entry digest, and precedence. Disabling retains the global block with
`enabled: false` and restores ordinary RMG rate sources and coefficients.

The fast numeric check independently reproduces these values with RMG's R:

| T (K) | PLP-SEC k (L/(mol*s)) | sensitivity RMG / PLP-SEC |
| --- | --- | --- |
| 600 | 63257.201921 | 0.130153877 |
| 700 | 160434.678657 | 0.183633454 |
| 800 | 322432.976207 | 0.248214469 |

The sensitivity check evaluates the previously recorded Arrhenius expression
A=926 cm^3/(mol*s), n=2.41, Ea=31547.4 J/mol on a mapped mock reaction.
It is not a new database estimate or an event-set compile.

## Phase-2b checks

The existing default `benzylicEndTest.py::test_compiled_benzylic_declarations_and_anchors_require_phase2b`
now also calls `plpsecLibraryTest.assert_compiled_plpsec_library`. The latter
independently enumerates the two expected complete graph pairs, requires
exactly one record per pair and **exactly two** direct library records overall,
and checks SI coefficients at **every** artifact grid point against the literal
Arrhenius expression. It also checks citation, global provenance, displaced
estimate tables, reciprocal links, forward/reverse units, and kf/kr=Kc.
Small inventories exercise the check; controls reject a doubled rate and an
extra library-tagged look-alike. It still fails on the old artifact at its
existing missing-declarations assertion; it cannot trigger a compile.

The existing phase-2b slow detailed-balance check
`compilerRealTest.py::test_all_pairs_have_exact_graphs_maps_degeneracies_and_detailed_balance`
now recognizes either an RMG family estimate or a library source as the directly
specified rate, then independently reconstructs thermo and verifies the inverse.
The manager's original post-compile default and focused slow commands thereby
cover this library as well. These future real-artifact checks were not run here.

## Reproduction

Run from `/home/alon/Code/RMG-Py-kmc-i050-plpsec-library`:

```bash
export PYTHONPATH=$PWD PYTHONHASHSEED=0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RMG_DATABASE_PATH=/home/alon/runs/i046-rules-from-training/database
export RMG_DATABASE_SHA=4a12d36fcdc193ede82c8d1ab5c1653495d445bc
export RMG_KMC_ARTIFACT=/home/alon/runs/i046-rules-from-training/cache/artifact/7491ed418f3ce2f0633688a9487cd12da278704b8f57e4f595176ddf80fb0296.json
export RMG_KMC_ALLOW_STALE_ARTIFACT=1 RMG_KMC_SLOW=0 RMG_KMC_PLPSEC_LIBRARY=1
export MPLCONFIGDIR=/tmp/i050-mpl COVERAGE_FILE=/tmp/i050.coverage
P=/home/alon/anaconda3/envs/rmg_env/bin/python
S=/home/alon/runs/i050-plpsec-library
ulimit -v 8388608
I050_CPUS=$("$P" -c 'import os; print(",".join(map(str, sorted(os.sched_getaffinity(0))[:4])))')
taskset -c "$I050_CPUS" "$P" -m pytest test/rmgpy/kmc/plpsecLibraryTest.py \
  test/rmgpy/kmc/benzylicEndTest.py test/rmgpy/kmc/compilerTest.py \
  test/rmgpy/kmc/cacheProvenanceTest.py test/rmgpy/kmc/rulesFromTrainingTest.py \
  -m 'not phase2b' -q -o addopts='' -o cache_dir=/tmp/i050-pytest-cache \
  --junitxml="$S/focused.xml" \
  > >(tee -a "$S/stdout.log") 2> >(tee -a "$S/stderr.log" >&2)
taskset -c "$I050_CPUS" "$P" -m pytest test/rmgpy/kmc -q -o addopts='' \
  -o cache_dir=/tmp/i050-pytest-cache --junitxml="$S/default.xml" \
  > >(tee -a "$S/stdout.log") 2> >(tee -a "$S/stderr.log" >&2)
```

The required `python setup.py build_ext --inplace` completed with exit 0
using the copied solver settings. The focused command above reproduced:

```text
87 passed, 1 deselected in 53.02s
```

This includes all 25 new library cases. The deselected test is the existing
phase-2b real-artifact acceptance check. GNU time recorded 54.79 s wall time
and 1,747,184 KiB peak RSS for this process. The explicit old artifact's
SHA-256 was independently reproduced as its filename, and the database HEAD
matched the dispatch's pin. AST parsing, `git diff --check`, and `bash -n`
on the reproduction block passed.

The complete default command reproduced exit 1 with:

```text
1 failed, 235 passed, 48 skipped, 3 warnings in 676.73s (0:11:16)
FAILED test/rmgpy/kmc/benzylicEndTest.py::test_compiled_benzylic_declarations_and_anchors_require_phase2b
AssertionError: phase-2b compile required: benzylic declarations absent
```

JUnit confirms 284 cases, one failure, zero errors, and 48 skips. This is the
only failure: the explicitly required old artifact has no new benzylic
declarations. The phase-2b check fails before attempting library acceptance;
it is neither skipped nor treated as passing. All 25 library cases pass in
both runs. The three warnings explicitly identify diagnostic stale-artifact
use. Existing slow/optional cases account for the skips under RMG_KMC_SLOW=0.
The 10,000-step default ledger check completed. GNU time recorded 11:22.44
wall time and 2,847,236 KiB peak RSS, below the 8 GiB cap. All processes used
at most four CPU affinity slots, with BLAS/OpenMP/MKL limited to one thread.

Production source did not change during either verifier run. Logs, JUnit XML,
exit statuses, and timing/RSS summaries are under
`/home/alon/runs/i050-plpsec-library/` (`stdout.log`, `stderr.log`,
`focused.xml`, `default.xml`, `focused.exit`, `default.exit`, and the two
`*.time.log` files).

Remaining work belongs to the manager: the single phase-2b event-set recompile
covering both the benzylic-end changes and this library, followed by the
post-compile default and focused slow acceptance checks. No new real-artifact
match count, grid coefficients, inverse rates, or ceiling are claimed here.
No real event-set compile, slow suite, tree regeneration, database write, or
prohibited dataset access is part of these commands.
