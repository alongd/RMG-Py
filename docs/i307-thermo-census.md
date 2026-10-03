# Thermo acquisition and consumer census

Resolved Species thermo is authorized only by `get_thermo_data()` and its whole-library comparison. Unsupported derivation paths refuse before changing data. A library representation may supply its processed NASA form and explicit energy and heat-capacity limits; copying attached values does not grant provenance.

The site inventory covers Python and Cython fields, constructors, thermo API calls, literal dynamic attribute accesses and writer templates in `rmgpy/`, `arkane/` and `scripts/`. Each record identifies the path, scope, expression and occurrence count. The static regression fails for an added, removed or changed site, including additions within an existing function. The tables group sites by their public boundary and give a behavioral regression for each verdict. Owner-free mathematical values and deferred storage are explicitly accounted for below.

Inventory: [excitedThermoSites.csv](../test/rmgpy/data/excitedThermoSites.csv). Machine-readable rows: [excitedThermoRows.json](../test/rmgpy/data/excitedThermoRows.json). Static tests: [excitedThermoCensusTest.py](../test/rmgpy/data/excitedThermoCensusTest.py).

## Acquisition and mutation

| Row | Routes | Verdict | Sites | Behavioral test |
|---|---|---|---:|---|
| A1 | Attached values and state cache | One structural library-source comparison; unknown model classes refuse by name | 21 | [test_nasa9_attachment_cannot_hide_between_five_cp_samples](../test/rmgpy/data/excitedThermoSourceTest.py) |
| A2 | Database library acquisition | Exact state and explicit limits only | 23 | [test_resolved_thermo_is_from_exact_library_entry](../test/rmgpy/data/excitedThermoTest.py) |
| A3 | GA, HBI, ML and surface estimation | Named refusal | 35 | [test_resolved_thermo_estimators_refuse_with_named_error](../test/rmgpy/data/excitedThermoTest.py) |
| A4 | Lower thermo ring/group/adsorption helpers | Named refusal, including atom-only ownership | 14 | [test_lower_acquisition_census_refuses](../test/rmgpy/data/excitedThermoSourceTest.py) |
| A5 | Cp0 and CpInf acquisition | Library limits or named refusal | 20 | [test_missing_library_limits_refuse_without_formula](../test/rmgpy/data/excitedThermoSourceTest.py) |
| A6 | Thermo representation processing and E0 | Validate library representation before processing | 36 | [test_attached_energy_and_limits_must_match_whole_library_source](../test/rmgpy/data/excitedThermoSourceTest.py) |
| A7 | Solvation and lower solute heuristics | Named refusal for derived thermo | 14 | [test_exact_solute_does_not_authorize_derived_solvation](../test/rmgpy/data/excitedThermoSourceTest.py) |
| A8 | Arkane input and thermo/statmech jobs | Named refusal | 73 | [test_jobs_and_mutators_census_refuses_before_use](../test/rmgpy/data/excitedThermoSourceTest.py) |
| A9 | Statmech database and mode fitting | Named refusal | 7 | [test_lower_acquisition_census_refuses](../test/rmgpy/data/excitedThermoSourceTest.py) |
| A10 | QM/QMTP dispatch and cache | Named refusal | 30 | [test_jobs_and_mutators_census_refuses_before_use](../test/rmgpy/data/excitedThermoSourceTest.py) |
| A11 | Isotope copying and entropy correction | Named refusal | 17 | [test_isotope_generation_and_direct_entropy_correction_refuse](../test/rmgpy/data/excitedThermoSourceTest.py) |
| A12 | Pressure-dependence mode/energy derivation | Named refusal | 36 | [test_network_census_refuses_even_with_cached_energy](../test/rmgpy/data/excitedThermoSourceTest.py) |
| A13 | Kinetics training/reconstruction and solute TS | Exact library acquisition; derived reconstruction/solute TS refuse | 66 | [test_kinetics_derivation_census_refuses](../test/rmgpy/data/excitedThermoSourceTest.py) |
| A14 | Sensitivity and uncertainty mutations | Named refusal | 28 | [test_jobs_and_mutators_census_refuses_before_use](../test/rmgpy/data/excitedThermoSourceTest.py) |
| A15 | Model admission and v=0 declaration | Getter validation, atomic declaration and label reservation | 21 | [test_admitted_automatic_label_preemption_is_atomic](../test/rmgpy/data/excitedThermoSourceTest.py) |

## Consumption and export

| Row | Routes | Verdict | Sites | Behavioral test |
|---|---|---|---:|---|
| C1 | Species scalar thermo readers | Checked getter | 12 | [test_statmech_cannot_supply_resolved_thermo_without_library](../test/rmgpy/data/excitedThermoReworkTest.py) |
| C2 | Species and Configuration conformer E0 | Refresh from checked getter on each read | 40 | [test_configuration_refreshes_energy_after_state_mutation](../test/rmgpy/data/excitedThermoSourceTest.py) |
| C3 | Species and configuration partition/density/mode readers | Named refusal | 37 | [test_statmech_consumers_census_refuses](../test/rmgpy/data/excitedThermoSourceTest.py) |
| C4 | Species to Cantera and Cantera model mutation | Checked getter before export | 11 | [test_cantera_refuses_unmatched_attached_state](../test/rmgpy/data/excitedThermoSourceTest.py) |
| C5 | Both Cantera YAML writers and bulk coverage reads | Checked getter before raw reads | 20 | [test_export_and_filter_census_rejects_altered_source](../test/rmgpy/data/excitedThermoSourceTest.py) |
| C6 | RMS YAML and RMS object conversion | Checked getter for RMS dictionaries; native Julia export of resolved owners refuses | 10 | [test_export_and_filter_census_rejects_altered_source](../test/rmgpy/data/excitedThermoSourceTest.py) |
| C7 | Chemkin thermo writer | Checked getter | 3 | [test_export_and_filter_census_rejects_altered_source](../test/rmgpy/data/excitedThermoSourceTest.py) |
| C8 | Thermo library export and reload | Checked library source for permitted exports; preserve v=0 and energy/limits; upstream Arkane explicit-state export refusal retained | 30 | [test_manifold_library_export_preserves_v0_identity](../test/rmgpy/data/excitedThermoSourceTest.py) |
| C9 | Arkane records, conformer output and plots | Named refusal; thermo-only writers use getter | 48 | [test_jobs_and_mutators_census_refuses_before_use](../test/rmgpy/data/excitedThermoSourceTest.py) |
| C10 | HTML thermo tables | Named upstream identity refusal before HTML rendering; permitted owners use checked thermo | 58 | [test_export_and_filter_census_rejects_altered_source](../test/rmgpy/data/excitedThermoSourceTest.py) |
| C11 | Filtering, merging, comparison and source extraction | Checked getter before numerical/raw attached reads | 63 | [test_export_and_filter_census_rejects_altered_source](../test/rmgpy/data/excitedThermoSourceTest.py) |
| C12 | Kinetics/TST energy consumers | Named refusal for state-blind statmech | 37 | [test_kinetics_derivation_census_refuses](../test/rmgpy/data/excitedThermoSourceTest.py) |
| C13 | Plasma thermo provenance | Checked whole library source before charge provenance | 7 | [test_plasma_refuses_altered_resolved_thermo](../test/rmgpy/solver/excitedThermoProvenanceTest.py) |
| C14 | Surface coverage thermo | Checked getter; derived coverage corrections refuse | 11 | [test_surface_census_refuses_derived_coverage](../test/rmgpy/data/excitedThermoSourceTest.py) |
| C15 | Network cached energies and pressure-dependence output | Named refusal; configuration-energy drawing uses checked E0 | 99 | [test_network_census_refuses_even_with_cached_energy](../test/rmgpy/data/excitedThermoSourceTest.py) |

## Sites without thermo authorization

| Row | Routes | Verdict | Sites | Behavioral test |
|---|---|---|---:|---|
| N1 | Owner-free thermo/statmech/kinetics models and ESS parsers | Mathematical values without a Species owner; owned acquisition/consumption is guarded | 231 | [test_owner_free_thermo_models_preserve_fields_in_round_trip](../test/rmgpy/data/excitedThermoSourceTest.py) |
| N2 | Deferred storage, copying, readiness and database services | Carries data/identity without authorizing it; all numerical consumers validate | 80 | [test_state_mutation_reloads_exact_state_thermo](../test/rmgpy/data/excitedThermoReworkTest.py) |

## Boundary notes

Atom-only ring heuristics check weak ownership registered by state setters, graph population, graph copies and ring extraction. Atoms retain their original fields. Shared atoms are checked against every live resolved owner; garbage collection releases registrations. The Molecule weakref slot is runtime metadata and is absent from explicit transport payloads.

Clearing a state header or manifold declaration refreshes formerly resolved caches through the getter before export or energy consumption. Species conformer E0 is refreshed from the checked source at each Configuration read. Pressure-dependent statmech and solver methods refuse resolved owners, including direct compiled solvers with cached arrays. Reference and BAC enthalpy caches recheck the current stored adjacency identity before use.

N1 contains model arithmetic, transition-state data, ESS parsing and external Cantera records. N2 contains constructor/copy/serialization data, readiness checks and database loading services. These carry values without authorizing their use: owner-aware estimation and numerical/export boundaries independently validate or refuse. The census is a source inventory and regression guard, not a proof against arbitrary reflection or caller edits to internal arrays.

A declared manifold exports the library lookup structure with vibrationallevel 0. ThermoData serialization retains E0, Cp0 and CpInf. Unsupported solvation processing refuses even when an exact-state solute entry exists. A refused solvent request leaves the prior gas source usable.
