# Qualified confined-anion wall interfaces and declared axial extension

The implementation supplies explicit `confinedAnion` and
`electropositiveBracket` closures, with `fullFrequency` production and
`radialOnly` sensitivity. The full-frequency multiplier is a **FINITE-CYLINDER
GEOMETRIC EXTENSION**. Kemaneci et al., arXiv:1612.07268, derive the factor
radially with axial loss neglected. They do not derive the axial multiplier.

The radial normalization integral and complete radial expression are retained
in `verify_mapping.py`. Profile normalization cancels in the relative factor.
`mappingTest.py` replaces the former absolute finite-cylinder equality probes
with the owner's infinite-column, alpha-zero, no-axial and synthetic no-radial
limit tests. The latter is a model-extension test, not source validation.

For the frozen 0.05 m radius, 0.30 m length deck, the axial EP share is
`f_z = 0.04525379687210524`. The actual difference between the geometry arms
relative to the EP total is `(1-h)*f_z`; it is not the axial share itself.
`verify_mapping.py --output <path>` records all exact frozen alpha points,
frequencies, differences and both ion-temperature conventions.

Gate C-radial compares the complete source expression after centre-to-volume
normalization. Its fixed numerical allowance covers the existing engine's
rounded 2.405 root. It is never used as the independent Gate C reference.
Gate C-geometry and full-profile C accept external reference callables and
frozen thresholds. **No production reference, epsilon or Report-10 budget is
supplied by this change. Finite-anion confinedAnion states refuse without them.**
The references in the new tests are explicitly test-only fixtures.

Build with `PATH=/home/alon/anaconda3/envs/rmg_env/bin:$PATH make build`.
Run Python from that environment and assert its `rmgpy.__file__` resolves
under the tested tree before any run. The regressions are
`test/rmgpy/solver/plasmaElectronegativeWallTest.py`,
`test/rmgpy/rmg/electronegativeWallInputTest.py`, and this directory's
`mappingTest.py`. Input-file semantics are documented in
`documentation/source/users/rmg/input.rst`.

`verify_bit_identity.py --output <dir>` dumps complete binary wall flux,
latched frequency, energy budget and trajectory for argon map mode and both
closures in each geometry arm with a zero-population anion. The old API
refuses even a zero core anion; its reference uses an inert zero neutral row
to hold the packed solver dimension fixed. No base guard is disabled. Both
inactive rows are checked to stay exactly zero; use `cmp` on every artifact.

The continuation report and actual logs are in
`/home/alon/runs/i313-en-wall/continuation/`. Real envelope qualification,
Report-10 sensitivity/materiality runs and merged-tip verification remain
outside this dispatch. The branch implementation does not close those
scientific gates.


The Gate-C adapter context follows I-314 §9: actual declared cylinder radius and
length, both-sign evaluated transport, accepted reaction propensities and charged
reactant classifications are available; arbitrary reference diagnostics survive
in the manifest. Eigenvalue-only declarations are labelled and never presented
as a declared chamber radius. The reference study's reported disagreement is an
owner decision; this implementation retains the ruled formula and fails Gate C
when a supplied reference disagrees beyond the frozen threshold. No real
reference was installed or run here, and no threshold/budget was chosen.

The numerical publication domain is defined once in
`rmgpy.solver.electronegative.EN_WALL_DOMAIN` and documented in the
[input-file guide](../../documentation/source/users/rmg/input.rst). It bounds
species concentrations, electron temperature, the axial/radial eigenvalue
ratio, EOS inventory volume and absolute cation EP wall frequency. Public
wall Jacobians and accepted-state monitors refuse outside-domain states by
variable name before publishing a result. Existing transport and A/B/C
qualification gates further restrict acceptance.

Every accepted finite-anion energy state checks its complete energy Jacobian,
including chemical electron dilution. There is no periodic-check gap. The
charged-wall operator rescales inventories, computes both geometry shares
from their own terms, and cancels reciprocal-inventory differences analytically.
Blanc mobility gradients sum per-bath differences so trace bath contributions
are retained when the mean mobility rounds to the dominant bath value.
Mass-action chemistry, the EOS and supported electron-temperature rate laws
have analytic energy-mode derivatives; no perturbation floor is used for
finite-anion chemistry. Energy diagnostics obtain the floating potential
without constructing a dense wall matrix. The zero-anion legacy solver policy
and its potential multiplication order remain unchanged.

`plasmaElectronegativeRework4Test.py` sweeps 12,960 in-domain states at the
corners and interior, every wall row and column, both closures and geometry
arms, both temperature modes, one/two-cation mixtures and flat/Blanc transport. Its independent
450-digit Decimal dual-number implementation differentiates the declared constitutive
law, with relative tolerance 2e-10 and no absolute tolerance. It also covers
quadratic and linear trace chemistry, tiny axial shares, complete dilution
guards, named outside-domain refusals and volume/frequency endpoints.

The rework-2 and rework-3 tests retain their numerical accuracy checks within
the publication domain. Former subnormal/extreme-inventory acceptance cases
now check named refusals; finite trace charged-pressure checks remain active
above the concentration floor. The numerical derivative references and the
source-mapping limits are not scientific qualification of the closure or
axial extension. The qualification fixtures echo production flux: they test
interfaces, refusal and recording only. No real full-profile reference,
physical error budget or production campaign is supplied here.
