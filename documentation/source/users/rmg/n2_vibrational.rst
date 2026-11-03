N2 vibrational ladder artifacts
===============================

The offline API ``rmgpy.tools.eedf.n2_vibrational.generate`` builds a sourced
N2(X,v) ladder, retaining v=0 through K (default 10, maximum 30). It returns
bytes for ``manifest.json``, a species/manifold ``input.py`` fragment and an
unchanged supplied ``thermo.py`` library snapshot. No reaction rate is fitted,
no thermo is estimated, and no electron cross-section array is copied into the
artifacts. This does not generate an EEDF table or execute a plasma reactor.

Local source manifest
---------------------

The source manifest is JSON. Paths must be absolute. Each pin has ``path``,
``sha256``, ``citation`` and ``provenance: "sourced"``. A ZIP pin may also name
an exact archive ``member``; the generator hashes the original archive snapshot
before reading that member. Files are read and hashed as one snapshot.
``estimated`` is reserved for a future estimator and refuses
in this implementation. The following is a template, not physical input::

    {
      "schema": "n2-vibrational-source-v1",
      "level_source": {
        "path": "/local/source/levels.json",
        "sha256": "<64 lowercase hexadecimal digits>",
        "citation": "<term energies, weights, isotope and convention source>",
        "provenance": "sourced"
      },
      "cross_sections": ["<pin objects for local LXCat Laporta files>"],
      "source_domain": {"initial_levels": [0, 1, 2], "final_v": 58},
      "source_assignments": {"0": "laporta_batch1", "1": "laporta_batch1"},
      "thermo_library": "<pin object for a fixed-v RMG thermo library>",
      "energy_tolerance_eV": "<explicit allowed source threshold/term mismatch>",
      "D6": {
        "set_by": "owner",
        "top_fraction_limit": null,
        "top_leak_share_limit": null,
        "K_step": 5
      }
    }

``levels.json`` declares ``molecule: "14N2"``,
``electronic_state: "X1Sigma_g+"``, ``energy_unit: "eV"`` and
``energy_zero: "v=0"``. Its ``levels`` array, ordered contiguously from zero,
contains ``v``, ``energy_eV``, positive ``statistical_weight`` and
``thermo_entry`` (the label in the pinned thermo library) for retained levels;
source-only levels above K may use ``null``. Energies must strictly increase,
starting at zero. Primary term-value/weight records must support this
explicit convention; ensemble vibrational thermo cannot substitute for fixed-v
library data. Every retained level's thermo entry is loaded and checked against
its exact molecular identity through the library lookup, with no estimation
fallback. Each entry must contain a NASA, Wilhoit, or ThermoData model; scalar
or other placeholder values refuse generation. The library snapshot is trusted
RMG Python input; loading it executes its database declarations, as normal
library loading does.

Coverage separates ``source_domain.initial_levels`` from ``source_domain.final_v``.
Every upward pair for each covered initial level through the inclusive final
domain is required, including transitions above K. ``source_assignments`` maps
each covered initial level (string key) to the ``source_id`` of its supplying
pin; records for that initial level in other files are ignored. Every initial
level 0..K must be covered, so a K requiring an uncovered initial level refuses
with that level named. Missing pairs, duplicates, out-of-domain pairs,
unsupported identities, malformed grids, hash mismatches and thresholds
inconsistent with declared term energies refuse before artifact emission.
Zero-threshold diagonal records are optional and audited. Source files remain
local and must never be copied into a repository.

D6 closure and coupling seam
----------------------------

The real top level has ``vibrationallevel K``. There is no separate lump species,
no sentinel, and no source-resolved tail population. Excitations from retained
initial levels to source final levels w>K use the original source cross section
and full electron threshold, redirecting their product to v=K. Their per-event
``top-leak:N2:vi-vw`` ledger records E_w-E_K, with positive molecular and negative
electron entries, and zero gas heating. Those entries are equal and opposite
views of the same leak, not two additional energy charges. Source-threshold
rounding differences are separately exposed as ``source_energy_residual_eV``;
the original threshold is never adjusted. For each channel::

    electron loss = retained energy gain + top leak + source energy residual

Source pairs with initial v=K and final w>K remain top-to-top energy-only
collisions with the full source threshold and ledger; they do not manufacture
a population above K. Records whose initial level exceeds K are explicitly
listed in the audit as ``unresolved_initial_above_K``; there is no resolved
target for their runtime flux. They are reclassified on regeneration at higher K.
Every supplied record is in the audit; every above-K excitation from a resolved
target is in the channel list and ledger.

Only pairs with both ends retained and distinct get a separate irreversible
superelastic partner. Both directions bind to the same ``source_pair_index``
and source block: forward to LoKI ``k_ine``, superelastic to ``k_sup``. Detailed
balance uses the supplied level statistical weights from that same source pair,
never thermal reversal. The source-pair index is a deterministic artifact index,
not a promise about the row order of a future solver collision set. A coupling
adapter must map the pinned file/block to its solver row. The artifact is a
versioned generation manifest, not a direct replacement for the qualified EEDF
channel-map schema. It declares hierarchical selector support and D6 execution
as runtime requirements and always reports ``runtime_ready: false`` here.

``assess_truncation(K, populations, forward_vibrational_power, top_leak_power,
policy)`` is pure. It normalizes nonnegative populations for v=0..K and compares
the top fraction and top-leak power divided by gross forward vibrational
electron power against the explicit owner limits. With zero forward and leak
power, the share is zero; inconsistent, negative or nonfinite input refuses.
Exceeding either limit requests regeneration at ``min(K + K_step, 30)``.
Exceeding at K=30 raises ``LadderError`` naming the failed limits. Equality
accepts. The coupling ticket owns applying the decision, regenerating and
restarting; artifacts alone do not implement the runtime check.

No physical default is supplied for either threshold. ``null`` placeholders
refuse both production generation and runtime evaluation. ``production=False``
(or CLI ``--draft``) allows a nonproduction artifact with placeholders, but
never bypasses source completeness or exact thermo validation. The
``production_manifest`` flag describes complete generation inputs and set
limits, not qualified plasma-runtime readiness.

Use from Python, or create a new output directory with::

    python -m rmgpy.tools.eedf.n2_vibrational /local/source.json /local/output --K 10

The CLI refuses to overwrite an existing output directory. Install the emitted
library in a caller-owned database as ``N2Vibrational`` before using the deck;
the generator never edits a database. The emitted deck is a loadable
species/manifold fragment; electron channel mappings require the coupling
adapter before a reactor can run. Retained v=0 has an explicit library identity,
while the unkeyed input ``N2`` uses ``vibrationalManifold(species='N2')`` for its
ground lookup. No Treanor or heavy-particle rate path is provided.

API
---

.. automodule:: rmgpy.tools.eedf.n2_vibrational
    :members: generate, read_cross_sections, require_complete_pairs, assess_truncation,
              LadderError, MissingPairsError

Resolved-population EEDF coupling
---------------------------------

``rmgpy.tools.eedf.vibrational.ladder_coupling`` adapts a generated manifest
into an optional ``vibrational_ladders`` entry in the EEDF generation spec.
Pass ``axis='Tv_N2'``, ``parent='N2(X)'`` and a list of collision bindings,
each with ``source_sha256``, ``source_block`` and the exact LoKI
``description``. Every audited source block must be bound, including diagonal
and unresolved-initial records. Ordering is immaterial: generator pair indices
are never interpreted as solver column indices. The adapter reads and verifies
the pinned spectroscopy snapshot and records all source levels, including the
unresolved levels whose populations are zero. The resulting metadata, source
pins, owner-set D6 policy and axis are included in the table fingerprint.
This adapter does not convert a raw LXCat file to a LoKI collision deck; the
caller must provide that deck and its independently validated channel map.

The selector and reactor implementation is independent of the molecular gas.
Selectors ``N2(X,v=1)``, ``O2(X,v=1)`` and ``NO(X,v=1)`` identify neutral
vibrational children of explicit electronic parents. Charged, wildcard and
rotational selectors still refuse. Electronic identities use the engine's
exact molecular ``electronic_state``, multiplicity and ``vibrational_level``;
labels do not select a state. Construct retained species from the manifest's
``states[*].adjacency_list``, including explicit v=0. An unkeyed vibrational
manifold ground species does not silently stand in for that resolved v=0.
Core species above K refuse rather than having their populations zeroed.
Only the N2 generator is supplied here; no O2 or NO ladder data is added.

The spec's ordinary ``axes`` mapping declares ``Tv_N2`` in kelvin. Each
retained or unresolved source level gets an explicit state-property declaration,
for example::

    N2(X) = 1
    N2(X,v=0) = boltzmann:Tv_N2:0
    N2(X,v=1) = boltzmann:Tv_N2:1
    N2(X,v=2) = boltzmann:Tv_N2:2

The last integer selects a sourced level, not a solver column. During table
generation these declarations evaluate a Boltzmann distribution normalized over
v=0..K, with zero population above K. Declare sourced level energies and
statistical weights as ordinary numeric state properties. Parent electronic
populations are fractions within the gas; vibrational populations are fractions
within their electronic parent. Channel-map ``target_fraction`` and
``product_fraction`` may name ``state:N2(X,v=0)`` and ``state:N2(X,v=1)``;
these resolve to gas fraction times electronic-parent fraction times vibrational
fraction. Numeric and coordinate-placeholder fractions retain their existing
behavior. The ordinary ``k_ine`` and ``k_sup`` columns are unchanged.

The reactor defines the v=0/v=1 moment as
``Tv = (E1-E0)/kB / log(n0*g1/(n1*g0))``, using sourced energies and weights.
Both populations must be positive and finite. Missing/ambiguous levels,
zero population, ``n1 > n0``, or a nonpositive logarithm refuse explicitly.
In particular, there is no zero-temperature clamp or infinite-temperature
approximation. These restrictions also apply to table-generation points.
The analytic derivatives are ``dTv/dn0 = -Tv/(n0*log_ratio)`` and
``dTv/dn1 = Tv/(n1*log_ratio)``. They enter the existing composition energy
chain rule used by the DAE. Tv uses the generic table interpolation and domain
checks; even solver trials cannot silently clamp this coordinate. Numerical
Jacobian population perturbations use one-sided differences toward the
interior near a Tv domain boundary; out-of-domain accepted states still refuse.

Forward and separate irreversible superelastic reactions must use exactly the
retained endpoint species. Above-K excitations redirect to K; top-to-top,
diagonal and unresolved-initial collisions are energy-only B columns. No
superelastic reaction is manufactured from an unresolved level. The runtime
checks every bound column's identity and threshold and requires all expected
inelastic/superelastic owners. Terminal qualification substitutes exact level
populations for the table's Boltzmann approximation and retains the existing
direct-solver comparison; a successful integration alone is not qualification.

Before terminal qualification, D6 evaluates the exact accepted populations and
gross forward power using the same table rates. Exceeding a limit raises
``D6RegenerationRequired`` with ``D6_REGENERATE`` and the requested higher K;
at K=30 it raises ``D6CeilingExceeded`` with ``D6_REFUSE_K30``. Unset limits
refuse. Regeneration/restart is an explicit caller action; production runs do
not regenerate tables automatically. Successful decisions are recorded under
``vibrational_D6`` in the terminal state.

The full collision threshold is already charged in the existing electron
energy row, including the superelastic gain. ``vibrational_energy`` computes
gross forward power from target populations times ``k_ine`` times the sourced
full threshold, including transitions above K. Retained energy, top-leak energy
and source-rounding residual partition that gross power. The separate
``discretization_power`` records the table's net channel power minus its net
flux times the sourced threshold, including interpolation discrepancies; it
never enters the D6 denominator. Its named
``top_leak_ledger`` is diagnostic and is
never added again to ``Q_inelastic`` or gas heating. D6 separately uses actual
accepted target populations. This implementation and its tests use synthetic
fixtures; they do not qualify a physical molecular arm. Tables without
``vibrational_ladders`` retain the existing argon path.

See :doc:`../../reference/vibrational_eedf` for the coupling API.
