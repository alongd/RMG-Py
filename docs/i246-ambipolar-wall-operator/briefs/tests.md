# Brief: pytest coverage for the charged-particle wall operator

## Intent

A new charged-particle wall boundary operator was added to `PlasmaReactor` in
`/home/alon/Code/RMG-Py-i246-ambipolar-wall-operator/rmgpy/solver/plasma.pyx`. It is currently
covered only by a standalone harness script. Turn that coverage into a real pytest file that runs
in CI, and prove every test in it can fail.

## Verifier

```bash
export PATH=/home/alon/anaconda3/envs/rmg_env/bin:$PATH
cd /home/alon/Code/RMG-Py-i246-ambipolar-wall-operator
pytest test/rmgpy/solver/plasmaWallTest.py -p no:cacheprovider --no-cov -q
```

Expected: exit 0, all tests passing, zero skips.

AND, separately and non-negotiably, a file
`/home/alon/Code/RMG-Py-i246-ambipolar-wall-operator/docs/i246-ambipolar-wall-operator/logs/red-states.log`
containing, for EVERY test function you write, a logged demonstration that it can fail:

1. the name of the test,
2. the exact one-line edit you made to break the thing it pins (show the diff hunk),
3. the captured pytest output showing that test FAILING,
4. confirmation the edit was reverted.

A test with no logged red state is not done. This project has three times in five days shipped a
check that could not fail, so this is the part that is actually being asked for.

<!-- paste below this line -->

## Where the work is

- Worktree: `/home/alon/Code/RMG-Py-i246-ambipolar-wall-operator` (a git worktree; run everything
  from there, do not `cd` to any other checkout).
- The implementation: `rmgpy/solver/plasma.pyx`. Read `_configure_wall`, `_resolve_wall_state`,
  `compute_nu_wall`, `_apply_wall_terms`, `_apply_wall_jacobian`, `residual`, `jacobian`, and
  `set_initial_conditions`.
- A WORKING reference harness that already exercises all of this and passes:
  `docs/i246-ambipolar-wall-operator/verify_operator.py`, with its shared model builder
  `docs/i246-ambipolar-wall-operator/argon_wall_model.py`. **Read both first.** Your job is
  largely to port what they check into pytest, not to invent new checks.
- The existing plasma test file, whose fixture style and conventions you must match:
  `test/rmgpy/solver/plasmaTest.py`.
- `pytest.ini` sets `python_files = *Test.py`, so the new file MUST be named
  `test/rmgpy/solver/plasmaWallTest.py`.

## Environment

- The default shell `pytest` is the WRONG interpreter and dies on a NumPy ABI error before
  collection. Always prefix: `export PATH=/home/alon/anaconda3/envs/rmg_env/bin:$PATH`.
- `rmgpy/solver/plasma.pyx` is cythonized. If you change it (you should only need to, temporarily,
  for red states), rebuild before measuring:
  `python setup.py build_ext --inplace -j 8`. A stale `.so` runs the old code and will tell you
  your change did nothing.
- `./rmgrc` already exists and points at `../RMG-database-plasma/input`. That database directory is
  **READ-ONLY for this work**. Do not write to it. If something you need requires changing it,
  stop and report that instead.
- Capture both streams for every command you run, into
  `docs/i246-ambipolar-wall-operator/logs/`:
  `<command> > >(tee -a stdout.log) 2> >(tee -a stderr.log >&2)`. A crash prints its traceback to
  stderr only, and without it the trace is lost.

## What the tests must cover

One test function each, named for what it pins:

1. `nu_wall = D_a/Lambda**2` matches an independent closed-form re-derivation to machine precision.
2. Halving `Lambda**2` exactly doubles `nu_wall`; and `nu_wall * Lambda**2` is constant across a
   range of chamber radii (so the diffusivity carries no geometry).
3. The wall loss is first order in the **electron** population: `-(wall flux)/N_e == nu_wall`
   exactly. Include the assertion that normalising to the neutral population instead gives a
   number orders of magnitude away -- that mis-normalisation is a real error made earlier in this
   project and the test should say so.
4. Net charge flux to the wall is exactly `0.0` on a charge-neutral state, AND is non-zero on a
   deliberately non-neutral one (both directions, in the same test or two).
5. `d(net charge)/dt` from the full residual -- chemistry plus wall -- is zero.
6. Heavy atoms leave the gas at exactly `(1-gamma)*nu_wall*N_ion`, checked at `gamma` = 1.0, 0.5
   and 0.0.
7. A reactor with no wall declared produces residual and Jacobian values **bit-for-bit identical**
   to one whose wall parameters are absent -- use `==`, not `approx`.
8. Below threshold (electron temperature 0.5 eV, no external source) the electron population
   decays toward extinction and no negative species amount is produced anywhere on the way.
9. The analytic Jacobian agrees with a central finite difference. **Do NOT use a single step
   size**: a central difference has truncation error ~h^2 and roundoff ~eps/h, and at h/|y| = 1e-7
   every configuration -- including the pre-existing no-wall code path -- disagrees at ~1e-3 for
   reasons that have nothing to do with the wall. Scan h over `(1e-8 ... 1e-2)` and assert the
   MINIMUM relative disagreement is below 1e-5. Include a no-wall control arm. `jacobian_scan()` in
   `verify_operator.py` already does exactly this; reuse its shape.
10. An ionisation degree above `max_ionisation_degree` raises `PlasmaStateError` naming the
    ionisation degree; one below it does not raise.
11. A non-charge-neutral initial composition with `quasineutral_electron=True` raises
    `PlasmaStateError` naming the net charge. (Without this guard the solver dies with an
    unreadable convergence error instead.)
12. Declaring a diffusion length without an ion mobility, or vice versa, raises.
13. `wall_recycling` outside [0,1], a non-positive diffusion length, a non-positive mobility, and a
    negative ionisation source each raise.
14. With `quasineutral_electron=True`, the electron row of the residual is the net charge and
    carries no `dydt` dependence; and the same reactor integrates successfully.

## Rules

- Follow PEP 8 for new code. Do not run a formatter over anything. Do not touch existing files
  other than to add the new test file -- if you believe a change to `plasma.pyx` is needed, stop
  and report it rather than making it.
- Every source file in this repository carries the MIT licence header; copy it from
  `test/rmgpy/solver/plasmaTest.py`.
- Do not use `--no-verify`, do not skip a Cython rebuild to make something pass, and do not weaken
  an assertion to make it green -- if a check fails, report the failure.
- Do not push, do not merge, do not open a pull request, do not commit. Leave the changes in the
  working tree and report what you changed.
- No test may contain a target electron density -- not as a constant, not as a fixture, not as a
  threshold. Assertions are stated in loss frequencies, charge, ratios and conservation, never in
  an expected absolute electron density.

## Report back

Three sections: **Completed**, **Verification** (the verifier command you actually ran and its
real output, plus the red-state log), **Remaining work**.
