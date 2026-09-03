# i197 — reverse collision-limit check passes a pressure where a potential belongs

## The defect

`Reaction.check_collision_limit_violation` reverse branch called
`reverse_kinetics.get_rate_coefficient(condition[0], condition[1])`. `condition[1]` is a **pressure in
Pa**. For the charge-transfer reverse kinetics that `generate_reverse_rate_coefficient` can return, the
second positional argument is instead an **electrode potential in V**, so a pressure of order
`1e5 Pa` was evaluated as `1e5 V`. It does not crash; it returns a number the collision-limit
comparison then trusts — and the number underflows to `0.0`, so the reverse check had silently become a
no-op for charge transfer.

## The correct evaluation point (round-11 rework)

The first fix passed `V = 0`. That is wrong whenever the reference potential is nonzero. The
Butler-Volmer activation shift is `Ea -= alpha*electrons*F*(V - V0)`; it vanishes at the **reference
potential `V0`**, not at `0 V`. The forward branch already evaluates charge transfer at `V0`
(`reaction.py:1093`), and `generate_reverse_rate_coefficient` copies a possibly-nonzero `V0` onto the
reverse object (`reaction.py:1441`, `:1465`, both `kr = …(…, V0=(V0,'V'))`). Surface reversible
potentials of ~0.4 V appear in existing tests; with `F = 96485` that is a large spurious `Ea` shift.
So the reverse charge-transfer rate must be evaluated at **its own `V0`**.

Selection is by the `V0` **property**, not an `isinstance` allowlist (the allowlist is a drift defect
this repo has paid for before). `getattr(reverse_kinetics, 'V0', None)` both selects the branch and
supplies the value: present on exactly the charge-transfer types, absent/None on the pressure-taking
ones (measured below). A future potential-taking type carries `V0` by construction, with no allowlist
to keep in sync.

## Probe answers (measurements)

**Q1 — the argument really is a potential, at source.**
- `ArrheniusChargeTransfer.get_rate_coefficient(self, double T, double V=0.0)` —
  `rmgpy/kinetics/arrhenius.pyx:3105`; docstring: "`V` is a potential, not a pressure."
- `SurfaceChargeTransfer.get_rate_coefficient(self, double T, double V=0.0)` —
  `rmgpy/kinetics/surface.pyx:928`.
- `ArrheniusChargeTransferBM`'s second positional arg is an enthalpy in J/mol
  (`arrhenius.pyx:3407`) — but it cannot reach this call site (Q3).

**Q2 — magnitude of the error.** `ArrheniusChargeTransfer` A=1e10, Ea=20 kJ/mol, V0=0, alpha=0.5,
electrons=-1, T=1000 K:
- `V = 0`   → `k = 902253917.4394315` (≈ 9.0e8)
- `V = 1e5` → `k = 0.0`
The pressure-as-potential drives the Butler-Volmer exponent so far the reverse rate **underflows to
exactly 0**, and a zero reverse rate never exceeds any collision limit — the reverse check became a
silent no-op for charge transfer, worse than a crash.

**Q3 — which types reach the line with a real object.** `generate_reverse_rate_coefficient`
(`reaction.py:1474`) returns a `SurfaceChargeTransfer` (`reaction.py:1445`) and an
`ArrheniusChargeTransfer` (`reaction.py:1469`) for those forwards. `ArrheniusChargeTransferBM` and
`Marcus` are `KineticsModel` subclasses outside the dispatch chain, hit
`UnsupportedReverseRateError`, and are skipped by the reverse branch — so only the two charge-transfer
types arrive here, exactly the two the `V0` probe selects.

**Q4 — any other instance of the same substitution.** No. Forward branch (`reaction.py:2026`) routes
through `self.get_rate_coefficient`, which dispatches each charge-transfer type to its own reference
potential before the generic `(T, P)` fallback. The `reverse_*_rate` helpers evaluate charge transfer
at `V0` (`reaction.py:1444`, `:1468`) and Chebyshev at `P` (`reaction.py:1542`). All correct.

**Discriminator verification** (`getattr(x, 'V0', None)` and `.value_si`):
```
ArrheniusChargeTransfer   V0=(0.4,'V')  value_si=0.4
SurfaceChargeTransfer     V0=(0.4,'V')  value_si=0.4
Arrhenius                 V0=None       value_si=None
SurfaceArrhenius          V0=None       value_si=None
StickingCoefficient       V0=None       value_si=None
Marcus                    V0 absent (getattr default) -> None
```

## Verifier evidence

- `TestReverseCollisionLimitPotential`: the two charge-transfer cases (nonzero V0=0.4) are RED before
  the rework — recorded `[0]` (the round-10 `0`), asserted `[0.4]` — and GREEN after; the non-CT case
  is a regression guard (never red), recording the pressure `[1.0e5]`.
- `test/rmgpy/reactionTest.py`: baseline **87 passed / 2 skipped**, after **90 passed / 2 skipped**
  (+3 new tests).
- 5 torr argon deck: `MODEL GENERATION COMPLETED`, core **3 species / 1 reaction**, `PYTHON_EXIT=0`
  captured from the interpreter. (The `ElectronCollisionPlasma` reverse-skip warning is the expected
  `UnsupportedReverseRateError` path, not this fix's branch.)
