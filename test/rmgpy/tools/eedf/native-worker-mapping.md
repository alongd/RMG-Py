# Native worker regression mapping

Frozen source: `03302e24f9ab5e43ef9c3da9821e0a5e36b6e472`. R11.9 guard inventory and R11.10 L3 authorize replacement of mechanism-only tests. Public signal, application-child, scientific, artifact and deadline obligations remain; no test is removed to satisfy a process cap. The 2026-10-24 plasma exec amendment specifically retains the inherited `[exec]` node.

The standalone `_native_lifecycle_probe.py` matrix (raising-handler, audit-logger, finalizer-lock, slow-unpickle, sigstop, sigkill) and `_native_attribution_probe.py` before/after run serially without an extra pytest parent. Their results are mandatory handback evidence. The two retained storm tests still run 20 fresh owners each; the final 1000-fork qualification is a separate opt-in gate.

| Frozen test | Disposition | Current obligation / coverage |
|---|---|---|
| `test_two_different_scan_paths_are_not_averaged` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_reactor_scan_finds_stable_and_unstable_states_and_fold` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_failed_seed_is_reported_not_dropped` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_conservation_modes_are_removed_before_stability_is_decided` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_nonconservation_zero_mode_is_marginal_not_stable` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_step_halving_disagreement_is_marginal_not_admissible` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_continuation_timeout_is_recorded` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_plasma_adapter_runs_engine_and_reads_engine_discharge_class` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_conservation_verification_uses_resolved_perturbation_ensemble` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_scan_refuses_more_than_six_worker_processes` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_continuation_follows_cubic_branch_and_stops_at_bracketed_fold` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_continuation_checks_first_interval_for_an_eigenvalue_crossing` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_continuation_scale_aware_tangent_brackets_rescaled_cubic_fold` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_continuation_certification_brackets_paired_stability_crossing` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_continuation_certification_finds_off_center_stability_excursion` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_default_certification_grid_matches_the_declared_minimum_step` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_fast_single_valued_branch_crossing_is_not_labeled_fold` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_steep_smooth_stable_branch_completes_without_a_false_fold` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_smooth_state_extremum_is_neither_a_fold_nor_unresolved` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_continuation_snaps_roundoff_sized_remainder_to_exact_target` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_continuation_refines_a_tangent_departure_without_an_eigenvalue_signal` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_continuation_refines_a_hop_seen_only_in_a_trace_component` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_eigenvalue_motion_bound_refines_toward_a_narrow_excursion` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_fast_companion_component_cannot_dilute_a_branch_jump` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_tangent_reversal_is_refined_and_reported_unresolved_not_fold` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_continuation_never_steps_below_the_declared_floor` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_corrector_noise_is_not_read_as_a_tangent_reversal` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_hopf_crossing_is_labeled_hopf_not_fold` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_worker_timeout_and_post_result_shutdown_are_deadline_bounded` | Remapped, retained node | Same integration-timeout bound. N1-b now rejects a sent result when bootstrap never completes normally; no artifact is produced. |
| `test_continuation_carries_the_full_state_between_corrector_steps` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_crashed_worker_is_recorded_and_sleeping_sibling_is_reaped` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_a_killed_worker_is_reaped_before_shutdown_returns` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_a_worker_that_outlives_sigkill_is_logged_not_waited_on` | Remapped, retained node | Bounded private reap now raises NativeWorkerStatusError. The retained exec test verifies the scan-level exact pending-PID log and artifact preservation. |
| `test_sigint_inside_worker_bookkeeping_still_reaps_every_started_worker` | Remapped, retained node | Same public obligation; private native context/record hooks replace removed fork/reaper internals. |
| `test_workers_run_with_the_callers_handlers_and_mask_and_the_scan_restores_them` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_a_worker_reaped_elsewhere_is_not_reported_alive` | Remapped, retained node | A caller's SIGCHLD reaper leaves multiprocessing believing the worker runs. |
| `test_sigint_delivered_through_another_thread_still_reaps_every_worker` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_default_sigterm_reaps_every_worker_then_terminates_with_its_status` | Retained; synchronized signal boundary | Same two signal cases and exact-reap/status -15 assertions; second signal is delivered at observed cleanup entry, not after a fixed delay that may outlive the owner. |
| `test_a_python_sigterm_handler_is_chained_after_every_worker_is_reaped` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_default_sigterm_during_worker_start_reaps_the_worker` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_sigint_entering_cleanup_after_an_exception_still_reaps_every_worker` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_cancelling_a_nested_scan_leaves_no_worker_at_either_level` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_workers_exit_when_their_owner_is_killed` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_a_worker_dies_with_its_owner_whatever_runs_in_it` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_a_worker_dies_with_an_owner_left_unreaped` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_a_worker_whose_owner_died_before_it_was_armed_exits_at_once` | Retired mechanism | Removed synthetic prctl/wrong-parent helper test. Real owner death before the earlier child hook returns is retained in test_a_worker_whose_owner_died_in_an_earlier_hook_ends_within_5_s_of_it and the six-case standalone matrix (sigkill). |
| `test_foreign_and_restart_dispositions_are_byte_identical_after_a_scan` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_a_default_sigint_set_behind_pythons_back_still_reaps_every_worker` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_rdkits_leftover_sigint_handler_is_left_to_act_as_it_would` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_a_start_interrupted_after_its_fork_still_reaps_the_child` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_the_fork_hook_never_closes_a_descriptor_it_did_not_open` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_a_native_sigint_handler_is_never_written_at_the_first_scan` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_python_handlers_installed_before_the_first_scan_are_deferred` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_without_a_symbol_table_python_handlers_are_left_alone` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_a_default_signal_teaches_the_python_handler_without_touching_a_native_one` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_the_symbol_table_locates_the_live_python_handler_and_only_it` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_handlers_that_disagree_are_not_taken_for_pythons` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_a_default_hidden_behind_a_cached_ignore_still_reaps_every_worker` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_rdkits_reset_behind_an_ignored_sigint_still_reaps_every_worker` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_an_interrupted_start_reaps_its_child_however_long_an_earlier_hook_runs` | Remapped, retained node | Native identity survives a real raising SIGUSR1 inside start, while an earlier application child hook sleeps, renames the worker, or closes descriptors. Ordinary children forked by Python/native code, other threads or handlers survive. The exact worker is reaped before scan raises without waiting for the hook. |
| `test_reclaim_kills_only_an_identified_worker_and_releases_the_rest` | Retired mechanism | Native PID publication and ordinary app survival: retained interrupted-start and reclaim-edge cases; standalone attribution-before/after. Nonces/candidate reading are removed. |
| `test_a_worker_that_dies_during_its_start_is_reaped_not_left_a_zombie` | Remapped, retained node | A worker exits in an early or late child hook. Its C-recorded PID remains authoritative before Python registration returns; exact cleanup leaves no worker zombie and preserves the unreadable application child. |
| `test_an_undecided_child_is_released_and_logged_at_the_reclaim_deadline` | Remapped, retained node | Native-plan remapping of the inherited identity edge cases. Exec replacement ends the watchdog guarantee: bounded typed refusal and exact pending-PID logging preserve prior artifacts while the replacement remains alive. Other workers have C-recorded identity and are reaped even when they exit early or share ancestry with unreadable application children. |
| `test_an_application_child_forked_in_the_window_survives_when_no_worker_is_forked` | Remapped, retained node | An application handler/thread forks before the worker fork is interrupted or fails. The native record has PID zero, so cleanup never signals or reaps the application child. Fork failure is a typed native startup refusal. |
| `test_a_worker_held_by_a_sleeping_hook_is_released_at_the_deadline_then_reaped` | Remapped, retained node | The native watchdog is already live when an earlier child hook sleeps. Native cancellation and exact reap complete without waiting for a report or for the hook to return. |
| `test_an_unreadable_held_worker_is_reaped_by_its_late_report` | Remapped, retained node | Native identity/reaping remains exact with unreadable child memory. |
| `test_a_later_interrupt_does_not_restart_the_cleanup_deadline` | Remapped, retained node | A real raising signal after C cleanup cannot create a second budget. |
| `test_an_interrupt_in_the_release_warning_preserves_it_without_extending_cleanup` | Remapped, retained node | Exec-limit logging retries within the already spent native deadline. |
| `test_release_resumes_after_raising_application_signals` | Retired mechanism | native interrupt modes cleanup-before, cleanup-after and bookkeeping; native record survives lost-popen. No release handoff or warning prerequisite. |
| `test_release_resumes_at_commit_and_acknowledgment_bytecodes` | Retired mechanism | native exact reap cache/lost-popen and cleanup-before/after/bookkeeping; C atomics replace the nine Python publication/acknowledgment stores. |
| `test_release_survives_a_signal_in_the_retry_bookkeeping` | Remapped, retained node | The C guard cancels and exact-reaps before public unwind when a third signal escapes Python retry bookkeeping. |
| `test_release_consumer_survives_a_signal_after_logging_lock_acquisition` | Retired mechanism | standalone finalizer-lock and audit-logger: native exit proceeds independently of all Python logging locks; caller resumes and exactly reaps. |
| `test_release_handoff_never_waits_on_an_application_before_fork_hook` | Retired mechanism | retained stalled-owner before-fork test with an already active native worker; native hook ordering and standalone finalizer-lock. No handoff lock remains. |
| `test_consumer_adopts_when_cleanup_never_stores_a_deadline` | Retired mechanism | native interrupt bookkeeping and expiry protocol: no cleanup deadline is required for the worker lease; owner resumes and exactly reaps. |
| `test_a_healthy_scan_renews_workers_across_several_run_leases` | Remapped, retained node | Same public obligation; private native context/record hooks replace removed fork/reaper internals. |
| `test_a_stalled_scan_owner_cannot_return_adopted_results` | Remapped, retained node | Same public obligation; private native context/record hooks replace removed fork/reaper internals. |
| `test_release_survives_twenty_fresh_signal_storms` | Remapped, retained node | Twenty fresh owners under real 2 ms raising SIGUSR1, native exit and exact reap. |
| `test_release_survives_twenty_storms_continuing_past_unwind` | Remapped, retained node | The C guard cancels and exact-reaps before public unwind; each real storm continues for 200 ms afterward without another reap. |
| `test_release_warning_is_written_when_the_owner_exits_immediately` | Remapped, retained node | Native cancellation/reap precedes owner exit; no release-warning guarantee remains. |
| `test_a_worker_can_log_while_the_prefork_consumer_holds_a_logging_lock` | Remapped, retained node | An ordinary parent thread holds a handler lock across native worker fork. |
| `test_interrupted_consumer_preflight_starts_no_worker_and_can_retry` | Retired mechanism | native missing-handshake/create-error plus artifact-startup: explicit admission refusal and exact reap, no fallback consumer. |
| `test_preflight_retries_a_consumer_retiring_before_lease_admission` | Retired mechanism | native per-worker readiness handshake and startup-failure refusal replace shared consumer preflight/retirement. |
| `test_one_blocked_emitter_cannot_block_two_reaps_or_spawn_another_emitter` | Retired mechanism | standalone audit-logger/finalizer-lock plus two-worker healthy/stalled controls. No native emitter exists and no logging precedes native reaping. |
| `test_the_deferred_reaper_never_signals_a_child_that_has_not_reported` | Retired mechanism | retained application-child-before/no-worker and unreadable-window controls, standalone exact attribution. Reports no longer grant authority. |
| `test_the_deferred_reaper_reaps_a_late_report_by_its_exact_pid_then_exits` | Retired mechanism | retained sleeping-hook and unreadable held-worker native tests. Native PID is available before hooks; exact reap needs no report. |
| `test_a_reaper_that_cannot_start_leaves_the_child_alone_and_logs_it` | Retired mechanism | native create-error and missing-handshake typed refusals, no artifact and exact reap. No Python reaper is started. |
| `test_a_live_daemon_reaper_leaves_nothing_at_interpreter_exit` | Retired mechanism | retained immediate-owner-exit test and standalone SIGKILL; native watchdog exits without Python finalization. No Python reaper/emitter shutdown exists. |
| `test_an_interrupt_while_reaping_the_reported_worker_still_reaps_it` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_a_worker_whose_owner_died_in_an_earlier_hook_ends_within_5_s_of_it` | Remapped, retained node | An earlier child hook sleeps for 12 seconds. Killing the real owner stops its worker within five seconds, before that hook returns; the controller verifies the exact reap. |
| `test_a_nested_scan_delivering_early_still_reaps_the_outer_workers_first` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_matching_branch_uses_artifact_tolerances_and_refuses_nonunique_matches` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_matching_branch_accepts_only_solved_powers_not_recorded_segments` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_matching_branch_uses_recorded_certification_points_on_curved_path` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_interpolation_estimates_are_recorded_as_uncertified_diagnostics` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_algebraic_electron_stability_is_invariant_to_charge_row_scale` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_branch_scanner_has_the_project_mit_header` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_r220_case_a_parallel_branch_hop_is_never_accepted` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_r220_case_b_earlier_fast_stretch_cannot_loosen_a_later_step` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_r220_case_c_noisy_initial_resolve_cannot_widen_the_noise_floor` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_r220_case_f_smooth_shift_cannot_carry_the_path_onto_the_upper_branch` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_r220_case_g_first_step_of_a_leg_is_certified_at_real_settings` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_certificate_sees_branches_that_exchange_places_within_one_step` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_aliasing_at_the_last_gap_is_exposed_by_step_halving_and_refused` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_marginal_sample_stops_the_path_with_an_event` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_continuation_refuses_a_non_stable_initial_state` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_off_grid_power_is_refused_as_not_certified_with_certified_powers` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_reference_power_is_validated_and_reconciled_with_the_paths_at_load` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_scan_records_certificate_resolution_coverage_and_consistency` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_fold_bracket_is_unresolved_when_a_non_fold_failure_interrupts_it` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_fast_mode_stiffening_far_from_zero_is_not_refused_at_real_settings` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_adapter_missing_data_is_refused_not_defaulted` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_scan_and_continuation_refuse_missing_reference_power_or_seed_state` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_predictor_density_is_extrapolated_in_log_space_not_substituted` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_timed_out_seed_is_rechecked_and_recorded_in_the_artifact` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_steps_never_exceed_max_step_and_a_step_plus_floor_tail_is_split` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_a_short_tail_is_split_rather_than_taken_above_max_step` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_the_first_step_of_a_leg_stays_within_max_step_at_real_settings` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_eigen_event_bracket_is_unresolved_when_a_probe_fails` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_per_component_departure_sees_a_jump_a_fast_companion_would_dilute` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_limit_3_interior_single_gap_alias_is_declared_not_caught` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_aliasing_radius_is_recorded_per_path_absolutely_and_in_tolerances` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_scan_without_half_step_path_is_unconfirmed_and_refused` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_scan_refuses_continuation_without_a_declared_resolution` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_limit_4_first_step_of_a_leg_has_no_departure_trigger` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_initial_state_on_another_branch_than_its_seed_is_refused` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_chain_sample_verdicts_and_code_hash_are_recorded` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_last_chain_sample_verdict_is_enforced_although_the_endpoint_replaces_it` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_coverage_spans_both_paths_once_and_matching_uses_half_step_states` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_record_continuation_refuses_paths_certified_at_other_tolerances` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_matching_tolerances_edited_after_the_scan_are_refused_not_applied` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_limit_6_time_budget_stop_on_the_real_adapter_narrows_a_fold_bracket` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_limit_2_states_differing_only_off_identity_are_one_branch` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_limit_5_unstable_window_narrower_than_the_resolution_can_be_missed` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_matching_refuses_between_solved_powers_and_never_assigns_b_to_a` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_marginal_final_chain_sample_from_another_basin_stops_the_step` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_artifact_from_the_producer_before_the_final_sample_gate_is_refused` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_fresh_paths_carry_the_current_certificate_on_every_accepted_state` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_marginal_reference_branch_is_reported_but_never_matched` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_stable_reference_branch_without_a_certified_path_is_never_matched` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_marginal_reference_verdict_withholds_a_branch_even_with_certified_paths` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_branch_fixtures_carry_what_the_current_writer_records` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |
| `test_user_guide_states_the_certificate_and_its_limits_verbatim` | Retained | Same test body and parametrization; backing native orchestration replaces the old guards. |

## Native replacement faults and retained obligations

The old consumer/nonce mutations M1–M17 are not rerun as mutations of deleted code.
The following table records their surviving behavioral obligations and the current
regressions. Executed mutant results and their source copies live with the N1 handback.

| Old fault | Current behavioral check |
|---|---|
| M1 zero-copy candidate | Real before-fork handler application child survives; no native child record is invented. |
| M2 reported PID kills another holder | Real after-fork handler application child survives; native C PID differs from it; native N7 app-report fault is rejected. |
| M3 unreadable identity | Retained unreadable interrupted-start variants and native unreadable boundary; no /proc identity inference. |
| M4 unarmed start waits | Before-fork attribution has PID zero and bounded cleanup; startup error/missing-handshake refusal. |
| M5 sole nonce holder identity | N7 genuine application PID report cannot replace the fork return; application child remains live. |
| M6 kill unreported expiry child | Watchdog acts only on its own process; native cancellation sends no PID signal; external reap refuses acceptance. |
| M7 restart cleanup deadline | Retained exec-retry test observes the same terminal-to-cleanup bound after a real later signal. |
| M8 exception drops active worker | Lost-Popen and cleanup-before/after/bookkeeping interrupts retain native records and exact reap. |
| M9 premature handoff marker | No queue handoff; prelinked C record survives a real handler before Python assignment. |
| M10 fork-registry lock | Blocking child-hook and before-fork stalled-owner cases, plus held parent logging-lock control. |
| M11 no consumer adoption | N1/N2 late/missing watchdog faults fail; GIL and early-hook expiry do not require a Python consumer. |
| M12 done before handoff | N4 before-reap and N5 before-finalizer faults fail; normal commit is immediately followed by native exit. |
| M13 consumer logs directly | Standalone audit-logger and finalizer-lock matrix: worker exits while owner logging is blocked. |
| M14 warning before reap | Native watchdog never logs; caller cleanup caches reap before any diagnostic; same blocked-lock matrix. |
| M15 deadline None waits forever | Blocking hook expires on the native running lease with no cleanup deadline stored. |
| M16 no renewal | Healthy worker and healthy two-worker scan span several short leases; N3 expiry cannot be revived. |
| M17 adopted scan returns | Public stalled-owner and final-retirement tests refuse results/artifacts; N6 scientific-124 fault fails. |
| C1 consumer-only immediate-exit warning | Superseded by accepted native diagnostic limits; immediate owner exit still leaves no running worker. Exec cleanup retains an exact pending-PID error/log and never claims reap. |

The SIGINT-then-SIGTERM test now synchronizes its second signal with entry to
cleanup. Native cleanup may complete before the old fixed 20 ms sender gap;
sending TERM to an already-exited owner did not exercise signal priority. The
same node still requires both workers reaped and owner status -15. The single
SIGTERM case is unchanged. The final-retirement regression additionally pauses
an owner after normal completion/exact reap but before retirement and requires
`ScanLeaseExpired`, with prior artifact bytes and mtime unchanged.

## Unconditional caller-thread cleanup guard

The worker loop moved from `_integrate_many` to `_integrate_owned`, enclosed by
the native `guard` call. The bookkeeping trace still interrupts the assignment
`interrupted = exc` after two earlier real handler raises. Its target code object
now follows that move. The interrupt probe retains the escaping exception and
requires an exact cached reap and an absent process before its verdict; its old
post-unwind private `Owner.cleanup()` call is removed. A zombie is a failure.

`_native_cleanup_matrix.py` discovers executable statements in the complete
cleanup block and injects 1, 2, 3 and 5 real nested handler raises at each site.
It retains owners and tracebacks, observes the public reap/report outcome, and
only then rescues any pending test exec process. The saved incoming-source run
distinguishes existing controls from reproduced failures. The guarded run covers
the moved cleanup block; pending diagnostics are separately tested at their new
reporting function. These assertions do not depend on a later scan or collection.

`_native_report_fork_probe.py` runs standalone (owner, worker, application child)
for both a reporting-callback fork and a reporter-exception-destructor fork.
Real reporting signals then exercise C retry. The child must reject ownership,
retain its unreaped copy (a foreign wait would cache ECHILD), and make zero further
report calls or non-self signals. Parent controls either consume the naturally
exited worker's status exactly once or emit exactly one pending-PID diagnostic.
Test rescue occurs only after the public verdict. Run each of `callback-report`,
`destructor-report`, `callback-reap`, `destructor-reap`; no pytest supervisor is
used, preserving the three-process bound.

The blocking-hook probe preserves each hook's readiness timestamp in a test-owned
pipe. The parent checks the first record after exact reap, retaining every prior
assertion. CPython swallowing a child-hook assertion and a later hook overwriting
the native timestamp can no longer hide a watchdog started after the first hook.
The late-watchdog mutant must fail in the parent, not merely print an ignored
exception in the child.

## Application forks inside launcher callbacks

Run `_native_fork_hook_probe.py` standalone in `before`, `parent`, `child`,
`audit` and `control` modes, without an extra pytest parent. Real callbacks fork
an application copy while public `scan()` is on the stack. The copy must leave
with ownership refusal before its alarm, without executing the adapter,
changing its native record/cache, signaling or reporting. The original parent
must complete with a normal exact-reaped worker and an artifact. A parent-visible
pipe records any foreign bootstrap even if the copy later exits normally.
The maximum topology is owner, worker and application copy (three processes).
