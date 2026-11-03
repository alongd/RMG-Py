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
``sha256``, ``citation`` and ``provenance: "sourced"``. Files are read and hashed
as one snapshot. ``estimated`` is reserved for a future estimator and refuses
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
``thermo_entry`` (the label in the pinned thermo library). Energies must strictly
increase, starting at zero. Primary term-value/weight records must support this
explicit convention; ensemble vibrational thermo cannot substitute for fixed-v
library data. Every retained level's thermo entry is loaded and checked against
its exact molecular identity through the library lookup, with no estimation
fallback. Each entry must contain a NASA, Wilhoit, or ThermoData model; scalar
or other placeholder values refuse generation. The library snapshot is trusted
RMG Python input; loading it executes its database declarations, as normal
library loading does.

Coverage requires every upward pair within the entire declared source domain,
including above K, not merely transitions currently active in the retained
ladder. Missing pairs, duplicates, out-of-domain pairs, unsupported identities,
malformed grids, hash mismatches and thresholds inconsistent with declared term
energies refuse before artifact emission. Zero-threshold diagonal records are
optional and audited. The supplied local batch must be completed by its data
owner; this tool does not download or invent missing pairs. Source files remain
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
