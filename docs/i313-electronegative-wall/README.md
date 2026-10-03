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

A/B/C refusals retain their gate names. Before either closure publishes
last-valid or gate diagnostics, the complete operator must be finite.
Scalar potential diagnostics add no extra dense matrix beyond this guard.

Every accepted finite-anion state checks its complete solver Jacobian, in
both prescribed-Te and energy modes, including chemical electron dilution. There is no periodic-check gap. The
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


Reaction energy and net electron dilution are combined before normalisation.
One incident electron is cancelled in each directed mass-action product;
intermediate product exponents are carried separately. The EN residual and
analytic energy Jacobian share this expression. Public ion-wall components
and frequencies synchronise Te from the supplied state, refresh the EOS
volume after that synchronisation and check the publication domain. Internal
Newton wall evaluations use their separate constitutive operators.

`plasmaElectronegativeRework5Test.py` independently varies volume, uses an
unequal two-cation split and extends its 800-digit Decimal comparison to the
full energy row with attachment chemistry. It checks all wall columns in
both closures, geometry arms, temperature modes and mobility conventions,
with relative tolerance 2e-10 and no absolute tolerance. The reverse-rate
regression retains a representable threshold-weighted derivative when the
unweighted rate per electron would overflow. The finite-rate overflow case
checks that prescribed-Te monitoring preserves last-valid instead of
publishing an infinite operator. Domain bounds remain unchanged.

Reaction-coordinate EOS derivatives sum the other pressure contributions
before division, preserving trace contributions when one reactant supplies
almost all pressure. The historical electropositive finite-difference trials
retain one legacy energy expression on both sides of the inactive anion face;
accepted finite-anion residuals and analytic Jacobians use the reduced form.

Prescribed-Te chemistry also combines direct and EOS derivatives using the
other-pressure numerator. Finite-anion forward, reverse and network-leak
residual products carry their exponents separately, and Gate B includes
electron multiplicity in the protected attachment product. Public anion
transport data and frequencies apply the same supplied-state preparation
as the cation APIs.

Electron-temperature derivatives shift each rate law's analytic temperature
power before accumulation, retaining constant k*Te and EOS cancellations.
Elastic partner and source-partition derivatives sum the remaining pressure
terms directly. The argon arithmetic and domain bounds are unchanged.

`plasmaElectronegativeRework6Test.py` adds focused 800-digit Decimal witnesses
and a 200-digit complete-operator sweep over zero-cation and minority-cation
faces, the 30/70 split, and the equal-population saturation corner. It varies
volume independently and exercises elastic loss, absorbed power, external
ionisation and a Te-dependent ionisation law in its energy configurations.
Every species row, energy row, Jacobian column and residual is checked.
The qualification fixtures remain test-only; these comparisons verify the
declared equations, without supplying a physical reference or error budget.


Rework 9 parameter bounds and rate evaluation
-------------------------------------------

`EN_WALL_DOMAIN` also bounds signed reaction energies to `[-1e8, 1e8] J/mol`
and every evaluated temperature exponent `n` to `[-50, 50]`, including elastic
collision fits, composite Arrhenius components and cached Te refresh laws. Configuration and the shared
acceptance gate check live quantities and packed energy arrays, so replacing a
law or mutating its exponent cannot bypass a refusal. Errors identify the reaction,
value and interval. Existing domain bounds retain their values.

The standalone two-temperature evaluator checks the actual `self.Te` before
calling the law, including when no state vector is supplied. Its temperature
interval is the same `[0.05, 50] eV` used by the reaction-reference API.

The chemical energy derivative multiplies the rate coefficient into each signed
energy term before the temperature slope. FMA product errors remain in the sum,
which avoids an overflowing energy/slope intermediate when a tiny rate makes the
final derivative finite. Constitutive arithmetic remains available for private
trial calculations; public acceptance refuses parameters outside the bounds.

Evaluated rate-cache checks use typed loops without NumPy temporary arrays.
Temperature validation and refresh use C scalar finiteness checks. They
still read every live value at each gate, including direct cache mutations. The
retained paired benchmark compares 120 alternating pairs after four warmups in
both batch and five-endpoint modes, asserts identical final state, flux, frequency
and budget, and requires each median tip/base ratio to be at most 1.05.
