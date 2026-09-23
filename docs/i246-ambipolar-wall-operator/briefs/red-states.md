# Brief: finish the red-state demonstrations for the plasma wall tests

## Intent

`test/rmgpy/solver/plasmaWallTest.py` has 16 test functions and they all pass. Five of them have
a logged demonstration that they can FAIL; eleven do not. A test with no demonstrated red state
is not evidence — this project has three times in five days shipped a check that could not fail,
which is the entire reason this task exists. Produce the missing eleven.

## Verifier

```bash
export PATH=/home/alon/anaconda3/envs/rmg_env/bin:$PATH
cd /home/alon/Code/RMG-Py-i246-ambipolar-wall-operator
pytest test/rmgpy/solver/plasmaWallTest.py -p no:cacheprovider --no-cov -q
```

Expected at the end: exit 0, `20 passed`, and `git diff --stat` empty for
`rmgpy/solver/plasma.pyx` and `rmgpy/rmg/input.py` — every breakage reverted.

AND `docs/i246-ambipolar-wall-operator/logs/red-states.log` must contain a complete four-part
entry for **all sixteen** test functions, in the same format the existing five already use.

<!-- paste below this line -->

## Where the work is

Worktree: `/home/alon/Code/RMG-Py-i246-ambipolar-wall-operator`. Run everything from there; it is
a git worktree, do not `cd` to any other checkout.

- Tests: `test/rmgpy/solver/plasmaWallTest.py`
- Implementation: `rmgpy/solver/plasma.pyx` and `rmgpy/rmg/input.py`
- Existing log, whose format you must match exactly: `docs/i246-ambipolar-wall-operator/logs/red-states.log`
  — read its five existing entries first and copy their structure.

Everything is committed. `git status --short` should show nothing but untracked dotfiles before
you start and after you finish.

## The eleven that still need a red state

```
test_net_charge_drift_zero_under_full_residual
test_heavy_atom_recycling_conservation
test_no_wall_residual_and_jacobian_bitwise_reproducible
test_below_threshold_decays_to_extinction_without_negative_amounts
test_analytic_jacobian_matches_finite_difference_scan          (parametrized x5)
test_jacobian_scan_negative_control_can_fail
test_ionisation_degree_ceiling_raises_naming_ionisation_degree
test_nonneutral_initial_state_raises_under_quasineutral_electron
test_diffusion_length_and_mobility_must_be_declared_together
test_wall_parameter_bounds_each_raise
test_quasineutral_electron_row_is_net_charge_with_no_dydt_dependence
```

## The cycle, per test

1. Make ONE minimal edit to the implementation that removes or inverts exactly the thing that
   test pins — not a broad breakage that would fail many tests at once. The point is to show that
   *this* test detects *that* defect.
2. Rebuild: `python setup.py build_ext --inplace -j 8`. `plasma.pyx` is cythonized, so without a
   rebuild the `.so` runs the old code and the test will pass, telling you nothing. This is the
   single most common way this task is done wrong.
3. Run only that test: `pytest test/rmgpy/solver/plasmaWallTest.py::<name> -p no:cacheprovider --no-cov -q`
4. Capture the failure output.
5. Revert with `git checkout -- rmgpy/solver/plasma.pyx` (or `rmgpy/rmg/input.py`). That is the
   exact restore because the file is committed — do not hand-undo.
   **Never run `git checkout -- .` or `git checkout -- docs/` or any path-less form**: there is
   untracked work in this tree that it would destroy.
6. Rebuild again, re-run the single test, confirm it passes.
7. Append the four-part entry to `red-states.log`: test name, the exact diff hunk of the
   breakage, the captured failing output, and confirmation of the revert.

Keep stdout and stderr for every command:
`<command> > >(tee -a docs/i246-ambipolar-wall-operator/logs/stdout.log) 2> >(tee -a docs/i246-ambipolar-wall-operator/logs/stderr.log >&2)`

## Two that need thought rather than a one-liner

`test_no_wall_residual_and_jacobian_bitwise_reproducible` pins that a reactor with no wall
declared is unchanged. A breakage must therefore be something that leaks into the WALL-LESS
path — for example, applying a wall term without checking `has_wall`, or making
`charge_row_scale` non-unity when `quasineutral_electron` is False. Breaking something that only
runs when a wall exists will not make this test fail, and if you find that you cannot make it
fail with any minimal edit, SAY SO in the log and in your report: that is a finding about the
test, not a failure of the task.

`test_jacobian_scan_negative_control_can_fail` is itself a negative control — it asserts that a
deliberately inconsistent Jacobian IS caught. Its red state is therefore the inverse: make the
Jacobian consistent again where the test expects inconsistency, and show the test fails. Think
about what that means before editing.

## Rules

- Do NOT weaken, delete, skip or `xfail` any test to make something pass. If a test cannot be
  made to fail, report that; do not quietly change it.
- Do NOT leave any implementation file modified at the end.
- Do NOT commit, push, merge, or open a pull request.
- Do not touch `/home/alon/Code/RMG-database-plasma` — it is read-only for this work.
- Do not add a target electron density anywhere.

## Report back

Three sections: **Completed**, **Verification** (the verifier commands actually run and their
real output, plus a count of complete entries in `red-states.log`), **Remaining work**.
