Offline LoKI-B EEDF tables
=========================

The offline generator and :class:`rmgpy.solver.eedf.EEDFTable` share a
content fingerprint. They provide electron rates, transport, mean energy,
normalized EEDF and channel powers together. This facility does not yet
choose a reactor operating point. A production ``PlasmaReactor`` using a
table must additionally pass the terminal qualification described below.

Generate a table with the RMG Python environment::

    python -m rmgpy.tools.eedf.generate spec.yml

Tables, manifests, exact setup files and independent held-out verdicts are
stored outside git under ``output_root/<artifact_sha256>/``. ``generation_spec.json``
records the selected and refined axes; use that spec to reconstruct the model
identity for loading::

    from rmgpy.tools.eedf.schema import load_spec
    from rmgpy.tools.eedf.generate import model_inputs
    from rmgpy.solver.eedf import EEDFTable

    spec = load_spec('/tables/<artifact_sha256>/generation_spec.json')
    # Obtain this content address from the approved generation result.
    expected_artifact_sha256 = '<artifact_sha256>'
    table = EEDFTable.load('/tables/<artifact_sha256>', model_inputs(spec),
                           artifact_sha256=expected_artifact_sha256)
    row = table.row(u, composition, branch_id='branch_0')

``u`` is ln(E/N / 1 Td). ``composition`` maps each named composition axis to
its current value. It may also supply current values for every screened
non-axis envelope. An omitted envelope coordinate asserts its manifest
reference value; reactor integration must supply actual evolving values.
Repository SHAs are provenance only. Every cross section and property file,
channel-map bytes and mapped reaction identity, physical state, solver
options and interpolant source identity are checked. No load regenerates
a table. Named exceptions identify changed inputs, an unaccepted table,
extrapolation, envelope exits, ambiguous branches and folded or poorly
conditioned mean-energy dependence on u.

Specification schema
--------------------

All fields below are mandatory; unspecified policy has no numerical default.
JSON specs are also accepted. In YAML, write exponential numbers such as
``1.0e-6`` as numeric YAML values.

* ``schema_version``: 1; ``arm``: gas/additive/feed identity; ``Tg_K``, ``P_Pa``.
* ``shared_objects``: mapping of absolute build-library paths to SHA-256.
  Sanitized ``ldd`` resolution must match these pins before every solve.
  Both ldd and LoKI receive an environment built from an empty dictionary:
  PATH=/usr/bin:/bin, C locale, UTC, the declared OMP thread count and
  OMP_DYNAMIC=FALSE. No invoking-shell loader, allocator or locale switches
  survive. This environment is part of the model fingerprint.
* ``binary`` and ``cmake_cache``: each ``{path, sha256}``; ``loki_commit`` and
  ``compiler`` record the pinned C++ solver and its build. The binary is checked
  at driver creation and before every execution. The build is never modified.
* ``input_files``: mapping from scratch-relative names to
  ``{path, sha256, kind}``, with kind ``cross_section`` or ``property``.
  Paths are relative to the spec unless absolute.
* ``channel_map``: ``{path, sha256}``. The external JSON list must exactly match
  the ordered LoKI collision descriptions. Each entry has ``description``,
  ``kind``, ``classification`` (A/B/C/D), ``threshold_eV``, ``target_fraction``,
  ``product_fraction``, ``sigma_max_m2``, ``flux_group`` and ``reaction``.
  ``reaction`` is null or ``{library, index, repr}``. Attachment entries also
  contain ``cross_section: {energy_eV, sigma_m2}`` for their energy moment.
  Every entry supplies its tabulated ``cross_section`` for independent EEDF
  moment qualification; elastic entries supply the physical ``mass_ratio``.
  Fractions can be coordinate placeholders such as ``"{x_O}"``. Flux groups
  declare the species production/loss scale used in absolute rate tests.
* ``axes``: ``u`` and proposed composition coordinates, each a strictly
  increasing node list. Composition axes are screened before selection.
* ``screen``: ``{rule: any_quantity_exceeds, threshold_fraction, candidates}``.
  Each candidate supplies ``{name, reference, values}``. The generator measures
  direct perturbations at every primary u node against H1–H4, records the worst
  quantity and keeps a significant axis. An omitted significant axis refuses
  generation. Dropped axes require a covering envelope; every nonzero-width
  envelope requires a recorded screen. Screening one coordinate at a time
  does not establish the absence of coupled effects; terminal re-solves are
  the responsibility of the later reactor-coupling ticket.
* ``envelopes``: mapping to ``{reference, min, max}``; frozen unsupported
  populations should have zero-width envelopes.
* ``working_conditions``, ``gas_properties``, ``state_properties``: explicit
  LoKI properties. Strings can use named coordinate placeholders.
* ``solver_options``: ``eedfType: boltzmann``,
  ``ionizationOperatorType: usingSDCS``, ``growthModelType: temporal``,
  ``includeEECollisions``, and ``numerics`` (``energyGrid``,
  ``maxPowerBalanceRelError``, ``nonLinearRoutines`` including
  ``maxEedfRelError``, algorithm and mixing parameter).
* ``tolerances``: G1, H1, H2, H3, H4 each have ``rtol`` and ``atol``;
  H4 also has ``total_rtol``. Additionally, ``fingerprint_rtol``,
  ``normalization_atol``, ``F0: {norm: weighted_L1, atol: ...}``,
  ``channel_power_sum: {rtol: ..., atol: ...}``, and ``u_condition`` containing
  ``min_abs_denergy_du`` and ``min_scaled_denergy_du`` are required.
* ``floors``: ``rate_absolute``, ``eedf_dynamic_range``,
  ``rate_flux_fraction``, ``absolute_flux_fraction``,
  ``relative_power_share``, ``absolute_power_share``. Below-floor rates are
  flagged and checked absolutely. Rate interpolation uses ln(1 + k/floor),
  preserving identically zero channels with continuous below-floor behavior.
* ``held_out``: ``{lhs_count, seed, cell_midpoints: true}``. Independent points
  include each axis-cell midpoint at the other axes' nodes, joint cell
  midpoints and seeded Latin-hypercube points. The study scheme uses 32 LHS
  points; a smaller count is useful for a real-solver integration fixture.
* ``refinement``: ``{max_rounds}``. Failing cells receive third-point nodes.
  All old independent points stay held out; new subcell midpoints join them.
* ``repositories``: engine and database SHAs, recorded without comparison.
* ``output_root``, ``scratch_root``, ``omp_threads``, ``max_parallel``
  (1–5), ``timeout_s``. Each subprocess writes separate stdout/stderr logs.

Interpolation and qualification
-------------------------------

Tensor PCHIP reduces composition dimensions first and u last. Rate logarithms,
mean energy and all transport and power quantities use the same interpolation
call. The mean-energy derivative is the derivative of that exact interpolant.
Zero or sign-changing mean-energy secants refuse the branch; a small derivative
refuses the requested row. Nonphysical transport and EEDFs also refuse.

The generator compares ascending, descending and independent cold scans and
preserves distinct observed paths as separate branch tensors. It never
averages different EEDFs. Branch count or correspondence changes across the
composition tensor refuse generation. Every held-out cell/quantity failure
is recorded, and an unsuccessful qualification returns CLI exit status 2.
The caller must supply the approved artifact byte hash; the loader hashes and
reads one opened snapshot. Acceptance is derived from recorded verdicts.
Every rate-changing floor and all qualification policies enter model identity.
Training, scan and held-out grids agree bit for bit, including stored cell faces.
Symbolic populations, even constant expressions, refuse unless resolved to
numeric stored values. HF, adaptive and nonuniform grids currently refuse.
Production loads require acceptance; ``require_accepted=False`` is available
only for the offline generator's validation of a provisional artifact.

At an accepted reactor terminal state, qualification resolves the effective
gas fractions and electronic-state populations, renders one LoKI setup, and
freezes that exact setup with its solver and input-file hashes. LoKI executes
the frozen bytes without rendering again. Its output must report the same
execution identity before the direct result can be compared with the table
row. A passing comparison is staged with the setup and both rows; export is
admitted only after an acceptance manifest binding those artifacts is
published atomically. A missing, stale, failed, or tampered manifest refuses
export, even if an older file contains a ``PASS`` label.

Qualification keeps a successful candidate unpublished while restoring its
owned SIGTERM and SIGINT dispositions. A cancellation received by either
deferrer prevents publication, including during restoration when another
signal may already terminate the process. Only a clean release publishes
PASS, under the same attempt lock, with an exclusive hard link of the owned
staged inode as the success boundary. A known prior record is retired first;
between retirement and publication the record path can be absent. Absence is
non-PASS: ``require_qualified`` refuses it with an EEDF qualification error.
After release, signals follow their original dispositions; foreign
handlers, hard process death, and power loss retain their existing limits.

A returning handler may re-enter the provider after signal release. A
terminal failure then revokes suspended publication, including its staged
file, and remains authoritative when the outer call resumes. Completion
returns the current run's finalised record. Export continues to see a
candidate until publication commits; it cannot observe a premature PASS
latch. Failure bookkeeping repeats interrupted cleanup before propagating
an interruption, preserving the first finalised failure and its run identity.

``require_qualified(operation)`` returns a revocable admission for the current
run generation. Failure finalised during validation raises. Failure after
validation can revoke an admission already on the return stack; receiving an
object is not permission to export. Qualification records, status strings,
record paths, table manifests and the reactor's public qualification fields
remain diagnostic snapshots with their existing types and contents. A retained
snapshot may still say ``PASS`` after failure and cannot authorize output.
The reactor keeps the capability separately as ``terminal_admission`` on both
the already-published and newly-qualified terminal paths. A returning signal
handler can revoke it at Python's ``RETURN_VALUE`` boundary:
``validate_terminal_state`` can still return normally with a revoked capability.
**A normal return carries no acceptance or export authority.** Rechecking just
before returning cannot remove this limit. Reactor initialization clears its
terminal admission; detached capabilities remain bound to their original
provider generation and are revoked by its failure or supersession.

EEDF completion logs describe the end of computation and missing output, and
make no accepted-run or file-validity claim. Physical termination flags,
profiles, energy/wall diagnostics, model snapshots and historical records are
not substitutes for the capability. Non-EEDF completion wording is unchanged.

Future terminal-qualified exporters must prepare their bytes and use
``admission.commit_bytes(destination, payload)``. This is the only supported
route for consuming terminal export authority: one destination per call,
immutable bytes, no callbacks. Boolean conversion refuses. The method
revalidates the current generation and committed acceptance bytes, stages a
private file, and publishes with an atomic exclusive hard link to the same
inode. The destination must be absent and the filesystem must support hard
links; existing files are never overwritten. The staging name is then removed.
Failure or a newer qualification revokes both the admission and any suspended
commit, including a commit interrupted immediately after publication.

Revocation first truncates the owned inode through its open descriptor and
writes a non-PASS marker at offset zero. It never reopens a replacement path
for writing. Cleanup atomically captures the public name in a private
directory, compares the captured inode with the retained descriptor, and
unlinks only inside that private namespace. A foreign inode is restored by a
no-replace hard link, or retained and reported by path if restoration cannot
complete. Between capture and restore, that foreign file is briefly absent
from its public name. Filesystems without hard links fall back to
retain-and-report for foreign captures; they cannot publish qualified output.
An occupied restore target is never overwritten.

Private staging and capture directories are named ``.eedf-<uuid>.private``
and have mode 0700. They live in the public directory's enclosing directory
when that is on the same filesystem; at a filesystem boundary they live in
the public directory itself. Open directory descriptors anchor cleanup.
Each empty directory is removed with non-recursive ``rmdir``. Successful
qualification and export leave no private directory behind. A directory
holding retained foreign data is reported by its exact path and is never
recursively removed. This namespace is process-owned; it is not a defense
against a hostile same-UID actor accessing the private directory.

Descriptor retention follows live authority and active calls. After a completed
qualification, a provider retains at most one current acceptance owner and four
descriptors: the record stream, public parent, private parent and staging
directory. Supersession or failure explicitly retires that owner. Completed
exports retain no descriptors; settled terminal failure leaves zero owners and
zero descriptors. Historical snapshots and externally retained retired objects
do not extend resource lifetime. A later invalidation of a released path uses
the same prior-record adoption checks and releases the adopted owner afterward.

Suspended publication and capture operations keep their resources until their
actual return or unwind, including settlement retries after a second signal.
With A additional file owners and C capture scopes in flight, the conservative
bound is ``4 + 4*A + C`` descriptors and ``1 + A`` owners; A and C return to zero
when the operations finish. The native lifecycle helper drops the final use
pin after the Python body's return event. It closes all descriptor slots and
removes matching registry membership without calling Python or releasing the
GIL. Python signal handlers therefore observe either live resources or a fully
released owner. Repeated retirement cannot close a reused descriptor number.

If namespace removal is forbidden but descriptor invalidation succeeds,
cleanup refuses and reports the remaining non-PASS inode. The explicit
filesystem limit is forbidden removal together with forbidden inode
invalidation: refusal reports the surviving file but cannot guarantee removal
of its bytes. The permission-denial and foreign-inode cases are regression
tested, including second real signals during cleanup. Hard process death and
power loss cannot run Python cleanup; this is not a durable multi-file
transaction.
No current production exporter consumes this admission; existing EEDF-mode
Chemkin and Cantera writer guards continue to refuse those exports. Diagnostic
manifests and profiles do not claim terminal qualification.

The pinned study binary has two explicit limitations: it recomputes a linear
initial solution at each scan node, so these ordered scans do not prove seeded
continuation, and temporal-growth iteration counts are not emitted. Counts
are stored as -1 (unknown). Actual seed-controlled continuation and resolved
symbolic state ladders require a solver output/interface extension. Reactor
stability and continuation in absorbed power belong to later tickets. Demo
tolerances are fixture proposals, not frozen production materiality gates.

The manifest and loaded table expose ``branch_certification``. The pinned
solver currently records ``uncertified: unseeded scans`` for each branch; scan
agreement supplies no proof of uniqueness or stability. Held-out checks include
the probability-weighted L1 EEDF distance, attachment energies, distribution
consistency of mean energy, rates and powers, and per-channel/group power sums.

Physical-file and population integrity
--------------------------------------

State-property assignment syntax follows the pinned legacy parser: whitespace
must surround ``=``. Other strings, including ``population=1``, name files and
must be in ``input_files``. Recursive legacy and JSON state-property files are
hashed and flattened to explicit numeric assignments. Ambiguous prefixes,
duplicate resolved selectors and unresolved functions refuse. ``LXCatFilesExtra`` and
``effectiveCrossSectionPopulations`` are checked for file closure and refused
because their collision/density mapping is not yet supported.

Numeric neutral electronic-state populations are currently supported.
Wildcard, charged-parent and vibrational/rotational hierarchy selectors refuse
until their full state-tree reduction can be represented. Map target and
superelastic product fractions must agree bit for bit with the gas fraction
times the explicit electronic-state population written to LoKI. Unspecified
electronic states have zero population. The loader independently checks every
stored training and held-out row against the resolved setup constants.

For each collision group, the allowed power-sum discrepancy is
``min(atol, absolute_power_share * abs(field)) + rtol * scale``, where
``scale`` is the minimum of the grouped power, the materiality threshold
``relative_power_share * abs(field)``, and the nonzero material channel powers
in that group. The threshold cap also catches a channel that becomes smaller
than the threshold because it was undercounted. Materiality uses ``relative_power_share``.
Thus an absolute floor cannot swallow a material channel at low field, and a
large group cannot hide an error in a smaller material channel. The numeric
budgets remain spec fields; the demonstration uses rtol=0.005, atol=1e-24,
relative_power_share=0.001 and absolute_power_share=0.0001.

Collision metadata and identifiers
----------------------------------

The declared ``sigma_max_m2`` must equal the maximum of the stored cross
section; floors are calculated from that section. Collision identities, kinds,
thresholds and all cross-section samples are checked against the hashed
classic LXCat PARAM/comment/data blocks that LoKI actually reads. Elastic
mass ratios and ionization OPB parameters must match hashed gas properties.
Reversible channels accept LoKI's ``<->`` output and require positive explicit
setup statistical weights matching the reverse-rate ratio. Unsupported
collision metadata formats refuse qualification.

Class A requires an explicit mapped reaction; B is energy-only and has no
mapped reaction. C is diagnostic and D is unsupported: tables containing
these classes may be inspected offline with ``require_accepted=False``, but
production loading refuses them. The table returns the complete collision
set, so it cannot silently exclude a C or D column. State chemistry eligibility
of a B channel still belongs to the frozen mechanism mapping.

Declared flux-group labels are retained as provenance. Qualification bounds
absolute rate errors by each channel's own flux, a conservative bound on a
species total; changing a group label cannot loosen a verdict. All numeric
tolerances and floor policies remain explicit spec fields.

HDF5 mapping, axis and branch keys use reversible percent encoding with a
stored encoding marker. Original identifiers, including slash, brackets,
quotes, plus, star and spaces, are returned unchanged; JSON preserves names
as literal strings. Legacy simple HDF5 keys remain readable subject to the
implementation fingerprint check. Solver population selectors still obey the
supported explicit neutral-state grammar; arbitrary stored identifiers do not
expand that grammar. ``omp_threads`` and ``max_parallel`` require positive
integers (boolean values refuse), with at most five parallel solves.

Artifact storage and reverse-rate qualification
-----------------------------------------------

Publication writes the HDF5 file, manifest, generation spec and validation
report in a sibling temporary directory before renaming it into place.
A directory without a manifest refuses to load. HDF5 storage is inspected
before any dataset or attribute values are read: only internal contiguous
or chunked numeric datasets, built-in filters, unique hard links, and
numeric/string attributes are supported. External raw storage, virtual
datasets, soft/external links, plugin filters and references refuse. The
writer checks the same restrictions before publishing.

Every state-property selector (population, energy and statistical weight) is
resolved as the pinned LoKI parser resolves it. Equivalent explicit neutral
selectors such as ``Ar(1S0)``, ``Ar(1S0,)`` and ``Ar(,1S0)`` become one canonical
state key; duplicate resolved selectors refuse. Unsupported root, charged,
wildcard and hierarchical overrides refuse for all three properties.

A defined reverse coefficient above its numerical floor uses the H3 relative
tolerance against its own magnitude at each held-out point, including when
the product population is zero. Below-floor coefficients retain the declared
absolute criterion. Population-weighted power and flux checks still use the
resolved physical populations.
