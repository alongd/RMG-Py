# Argon resonance trapping implementation evidence

Implemented the explicit `holsteinWalsh` and `complete` reactor operators for
resolved Ar 3P1/1P1. See the user documentation for equations, geometry, the
refusal domain and the quantified approximation for omitted N2 broadening.

## Reproduced reacting evidence

The manager authorized read-only database review commit
`84b9596b3a53084564e0f2f70f7caaac463ec464`, which includes the higher-cascade
convention. Runtime build base: RMG-Py
`b634d0842593ffeb7b7d7d7bf197779414968fa3`, plus this feature.
The new vacuum-line JSON is from this feature's isolated database worktree.
No electron channel, cascade yield, metastable, N2/NO/O2 or manager-review-tree
file was modified.

A fresh direct LoKI-B solve used the production Ar converter, real cross sections
and all 39 collision descriptions. The conversion had zero missing/shared owners;
312 point consistency checks passed. All 33 cascade families had four owners and
one threshold energy debit per collision. The source binary, input and channel-map
hashes are recorded with the evidence.

For the reacting test, all 136 direct and cascade excitation owners were replayed
at that **fixed 60 Td EEDF**, with prescribed electron density 1e13 m^-3. These
are process-local copies with constant rate coefficients from the fresh solve;
the database EEDFChannel objects and yields are unchanged. Ionization and wall
loss are absent from this reduced subsystem. Four resolved 4s populations and
ground Ar evolve in the compiled PlasmaReactor. The two metastables continue to
accumulate; only the resonance subsystem is assessed for steady state. This is
not a table-interpolation qualification or a global-discharge sustainment result.

Conditions: 298.15 K, 666.61184 Pa (5 torr), radius 0.05 m, length 0.30 m,
EOS volume 0.002356194490192346 m^3. Both models reached 1 ms. The Holstein–Walsh
arm's resonance densities varied by less than 0.1% between 0.5, 0.75 and 1 ms;
feed/loss mismatch was below 0.1% at all three accepted samples.

| Model/state | Density at 1 ms (m^-3) | g | A g (s^-1) | Feed/loss relative mismatch at 1 ms |
| --- | ---: | ---: | ---: | ---: |
| Holstein–Walsh 1s4 | 2.666860241915321e15 | 2.887440031283162e-4 | 38114.20841294 | 3.33e-8 |
| Holstein–Walsh 1s2 | 5.689831308670779e14 | 2.870655513888687e-4 | 152718.8733389 | 8.34e-9 |
| Complete 1s4 | 1.016452677748976e17 | 0 | 0 | Not steady; zero radiative loss |
| Complete 1s2 | 8.689446546558194e16 | 0 | 0 | Not steady; zero radiative loss |

Holstein–Walsh photon power was 0.8340108789 W in that EOS volume; complete
trapping gave exactly zero. Both runs wrote ordinary Chemkin (136 excitation
replay reactions) and native CSV/radiation JSON. The operator itself was omitted
from Chemkin with the shared header annotation. Initial standalone writer attempts
required loading electron thermo and applying the normal NASA conversion; the
final recorded run completed both exports. Thermo conversion emitted SciPy
ill-conditioning warnings; this test does not use fitted thermo to derive the
excitation coefficients or radiation closure.

The database already cites Bogaerts, Gijbels and Vlcek (1998),
[doi:10.1063/1.368009](https://doi.org/10.1063/1.368009). Figure 2, p.130,
shows spatially varying resonance populations on a 10^12 cm^-3 scale at 1 torr,
1000 V and 2 mA, falling to zero at boundaries. The present conditional densities
(2.67e9 and 5.69e8 cm^-3) are much lower than its near-cathode populations.
The prescribed low electron inventory, absence of other excitation/loss pathways,
and different discharge conditions prevent interpreting that difference as a
validation error or treating the published profile as a 5-torr acceptance band.

## Reproduction artifacts

All commands used the dispatch's systemd 16 GiB/no-swap/800% CPU scope and
choom/nice/ionice. Worktree-local builds were used without changing the shared
Python installation. Full paper PDFs and cross-section payloads remain outside git.

Scratch root: `/home/alon/runs/argon-trapped-decay/`.

- `argon_point.py`: manager conversion/direct-solve recipe, with own output paths
  and 298.15 K; uses review database read-only.
- `real_reactor.py`: reduced reacting replay and export checks.
- `real-eedf/spec.json`, `conversion-report.json`, `channel-map.json`, `row.json`,
  `verdicts.json`: pinned setup and direct-solve records.
- `real-eedf/reactor-evidence.json`: feeds, source hashes and accepted balances.
- `real-eedf/holsteinWalsh/` and `real-eedf/complete/`: Chemkin and native profiles.
- `evidence/real-point.{stdout,stderr}` and `evidence/real-reactor.{stdout,stderr}`.
- `source-and-reactor-evidence.py`: checks the new vacuum-line data against the
  archived ASD response and exercises a separately labeled synthetic benchmark.
- `regenerate-censuses.py` and `evidence/census-regenerate.{stdout,stderr}`:
  explicit changed-site classifications, checked by the repository collectors.

Run the two evidence scripts with rmg_env Python and this Py worktree first on
PYTHONPATH, from this worktree (whose rmgrc points to its paired DB). The manager
reruns on the landed database tip at intake.

## Pre-review

- **Fixed:** an unresolved triplet Ar metastable was initially counted as an
  unsupported foreign gas. It now counts in Ar inventory while receiving no
  radiative loss; a compiled-reactor regression covers it.
- **Fixed:** exposed declaration dictionaries could be mutated independently of
  the live operator. Defensive property copies keep diagnostics, pickle and live
  model consistent; the mutation test failed before this fix and now passes.
- **Checked:** residual/Jacobian stoichiometry, ground-density/volume derivatives,
  exact complete endpoint, no second electron cooling, reinitialization/index
  rebinding, trial-state publication, physical geometry round-trip, source rows,
  actual output files, no vacuous steady-state assertions, and shared geometry
  validation rather than a second wall parser.
- **Rebutted:** nonfinite trial-state propagation is already refused by the
  existing EOS gate; the regression verifies that accepted diagnostics survive.
- **Not applicable:** no lock, concurrent publisher, retry/supersession or external
  transaction was introduced. Existing provider ownership was not changed.

Final worker tests: **289 passed, 1 skipped** (7.31 s), covering radiation,
compiled reactor integration, cascade fixtures, input parsing and all diagnostic
writers. Both build and source-row checks passed. `git diff --check` passed.
A separate RST parser check could not run because rmg_env lacks `docutils`;
no documentation build is claimed.

The full suites and full state-data/reversal/excited-export census tests are
manager-run. Worker evidence claims only the targeted tests and collectors.
