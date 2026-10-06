"""Generate the real table, prepare the pure-Ar deck, and run it when admissible."""

import argparse
import hashlib
import json
import logging
import math
import os
from pathlib import Path
import subprocess
import sys
import traceback

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from rmgpy import settings
from rmgpy.exceptions import PlasmaStateError
from rmgpy.rmg.main import RMG
from rmgpy.solver.eedf_provider import (
    DEVELOPMENT_UNQUALIFIED_STATUS,
    development_unqualified_table_route,
)
from rmgpy.solver.electronegative import manifest_values
from rmgpy.tools.eedf.loki import LoKIDriver, enrich_row
from rmgpy.tools.eedf.schema import file_hash, load_spec
from rmgpy.tools.eedf.validation import compare_row
from real_eedf_harness import (
    DEVELOPMENT_SWEEP_TOLERANCES,
    generate_real_artifact,
    development_argon_model,
    power_consistent_initialization,
    pure_argon_deck,
    reactor_readiness,
)


def development_diagnostic(blockers):
    """Return the mandatory warning fields for every harness diagnostic."""
    return {
        'scientific_status': DEVELOPMENT_UNQUALIFIED_STATUS,
        'artifact_blockers': list(blockers),
        'export_allowed': False,
        'qualification_allowed': False,
    }


def candidate_endpoint_classification(reactor):
    """Return the candidate label only for a measured steady-state stop."""
    if getattr(reactor, 'steady_state_reached', False) is not True:
        raise RuntimeError(
            'development run did not reach steady state; candidate endpoint refused')
    return 'candidate numerical steady state'


def write_json(path, value):
    """Write one strict, stable development diagnostic JSON file."""
    Path(path).write_text(json.dumps(
        manifest_values(value), indent=2, sort_keys=True, allow_nan=False))


def git_fingerprint(path):
    """Fingerprint a Git tree by commit plus the current tracked diff bytes."""
    path = Path(path).resolve()
    try:
        commit = subprocess.check_output(
            ['git', '-C', str(path), 'rev-parse', 'HEAD'], text=True).strip()
        diff = subprocess.check_output(
            ['git', '-C', str(path), 'diff', '--binary', 'HEAD'])
    except (OSError, subprocess.CalledProcessError) as error:
        return {'status': 'NOT EXECUTED', 'reason': str(error), 'path': str(path)}
    return {
        'path': str(path),
        'commit': commit,
        'tracked_diff_sha256': hashlib.sha256(diff).hexdigest(),
        'clean': not diff,
    }


def configuration_fingerprint(args, deck_path):
    """Fingerprint the deck and the runner settings that select this trajectory."""
    values = {
        'deck_sha256': file_hash(deck_path),
        'development_unqualified': bool(args.development_unqualified),
        'seed_multiplier': args.seed_multiplier,
        'arm_label': args.arm_label,
        'over_ionised_domain_guard': bool(args.over_ionised_domain_guard),
        'restart_from': args.restart_from,
        'progress_interval_seconds': args.progress_interval_seconds,
        'wall_budget_seconds': args.wall_budget_seconds,
        'sweep_tolerances': DEVELOPMENT_SWEEP_TOLERANCES,
    }
    encoded = json.dumps(values, sort_keys=True, separators=(',', ':')).encode()
    values['sha256'] = hashlib.sha256(encoded).hexdigest()
    return values


def select_restart_checkpoint(trajectory):
    """Select a measured intermediate accepted state near half convergence time."""
    if len(trajectory) < 3:
        raise RuntimeError('trajectory has no interior accepted state for restart')
    target = 0.5 * trajectory[-1]['time_s']
    interior = trajectory[1:-1]
    checkpoint = min(interior, key=lambda row: abs(row['time_s'] - target))
    result = dict(checkpoint)
    result['checkpoint_role'] = 'intermediate accepted state for restart test'
    return result


def a_posteriori_resolve(reactor, artifact, run_directory):
    """Directly re-solve LoKI-B once at the candidate endpoint and compare."""
    endpoint = reactor.development_trajectory[-1]
    result = {
        'check': 'a-posteriori LoKI-B endpoint re-solve',
        'artifact_sha256': reactor.electron_kinetics['table'][1],
        'candidate_state': endpoint,
    }
    spec_path = Path(artifact) / 'generation_spec.json'
    if not spec_path.is_file():
        result.update(outcome='NOT EXECUTED',
                      reason='artifact has no generation_spec.json')
        return result
    try:
        spec = load_spec(spec_path)
        spec['scratch_root'] = str(Path(run_directory) / 'loki-a-posteriori')
        channel_map = json.loads(Path(spec['channel_map']['path']).read_text())
        resolved_coordinates = reactor._eedf_coordinates(reactor.y)
        coordinate_names = (
            set(spec['axes']) | set(spec.get('envelopes', {}))) - {'u'}
        coordinates = {
            name: resolved_coordinates[name] for name in coordinate_names}
        driver = LoKIDriver(spec)
        direct = driver.run(
            coordinates, [float(reactor.eedf_row.EN_Td)], 'endpoint')[0]
        enrich_row(direct, channel_map, spec)
        direct['power_absolute'] = (
            abs(direct['power_groups']['field']) *
            spec['floors']['absolute_power_share'])
        predicted = reactor.eedf_row.as_dict()
        checks = compare_row(predicted, direct, reactor.eedf_provider.runtime_metadata)
        result.update(
            outcome=('PASS' if all(check['passed'] for check in checks) else 'FAIL'),
            checks=checks,
            coordinates=coordinates,
            table_prediction=predicted,
            direct_result=direct,
            comparison_tolerances={
                'tolerances': reactor.eedf_provider.runtime_metadata['tolerances'],
                'floors': reactor.eedf_provider.runtime_metadata['floors'],
            },
            direct_EN_Td=float(direct['EN_Td']),
            direct_power_groups=direct['power_groups'],
        )
    except (FileNotFoundError, KeyError) as error:
        result.update(outcome='NOT EXECUTED', reason=str(error))
    except Exception as error:
        result.update(outcome='FAIL', reason=type(error).__name__ + ': ' + str(error))
    return manifest_values(result)


def negative_fixture_result(reactor, error):
    """Classify the named old-seed domain-guard regression without softening it."""
    state = np.asarray(getattr(reactor, 'y', reactor.y0), dtype=float)
    message = str(error)
    checks = {
        'plasma_state_error': isinstance(error, PlasmaStateError),
        'named_table_domain_error': 'accepted EEDF state is outside the table domain' in message,
        'state_and_bound_identified': ('outside' in message and '[' in message and ']' in message),
        'last_accepted_state_finite': bool(np.all(np.isfinite(state))),
        'last_accepted_state_nonnegative': bool(np.all(state >= 0.0)),
        'provider_remained_loki_table': reactor.electron_kinetics['provider'] == 'loki-table',
    }
    return {
        'fixture': 'OVER-IONISED INITIAL-STATE DOMAIN-GUARD TEST',
        'outcome': 'PASS' if all(checks.values()) else 'FAIL',
        'checks': checks,
        'exception_type': type(error).__name__,
        'exception': message,
        'artifact_sha256': reactor.electron_kinetics['table'][1],
        'scientific_status': DEVELOPMENT_UNQUALIFIED_STATUS,
        'export_allowed': False,
        'qualification_allowed': False,
    }


def development_failure(reactor, blockers, error):
    """Describe a development-run stop, including the forced budget progress line."""
    failure = development_diagnostic(blockers)
    failure.update({
        'exception_type': type(error).__name__,
        'exception': str(error),
        'traceback': traceback.format_exc(),
    })
    if getattr(reactor, 'development_wall_budget_hit', False):
        failure['stop_reason'] = 'wall_budget'
        failure['last_progress'] = dict(reactor.development_last_progress)
    for output_name, attribute_name in (
            ('energy_budget', 'energy_budget'),
            ('latest_energy_terms', 'electron_energy_terms')):
        values = getattr(reactor, attribute_name, None)
        if values:
            failure[output_name] = manifest_values(dict(values))
    try:
        state = getattr(reactor, 'y', None)
        if state is None or not len(state):
            state = getattr(reactor, 'y0', None)
        if state is not None and len(state):
            coordinates = reactor._eedf_coordinates(state)
            failure.update({
                'final_time_s': float(getattr(reactor, 't', 0.)),
                'u': float(state[reactor.te_index]),
                'coordinates': coordinates,
                'Ar4s_fraction': float(coordinates.get('Ar4s_total', 0.)),
            })
    except Exception as diagnostic_error:
        failure['state_diagnostic_error'] = (
            type(diagnostic_error).__name__ + ': ' + str(diagnostic_error))
    return failure


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-directory', required=True)
    parser.add_argument('--artifact')
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--development-unqualified', action='store_true')
    parser.add_argument('--progress-interval-seconds', type=float, default=30.)
    parser.add_argument('--wall-budget-seconds', type=float)
    parser.add_argument('--seed-multiplier', type=float, default=1.0)
    parser.add_argument('--arm-label')
    parser.add_argument('--restart-from')
    parser.add_argument('--over-ionised-domain-guard', action='store_true')
    args = parser.parse_args()
    if (not math.isfinite(args.progress_interval_seconds) or
            args.progress_interval_seconds <= 0.):
        parser.error('--progress-interval-seconds must be finite and positive')
    if args.wall_budget_seconds is not None:
        if not math.isfinite(args.wall_budget_seconds) or args.wall_budget_seconds <= 0.:
            parser.error('--wall-budget-seconds must be finite and positive')
        if not args.development_unqualified:
            parser.error('--wall-budget-seconds is development-route only')
    if not math.isfinite(args.seed_multiplier) or args.seed_multiplier <= 0.:
        parser.error('--seed-multiplier must be finite and positive')
    if args.restart_from and args.over_ionised_domain_guard:
        parser.error('--restart-from and --over-ionised-domain-guard are mutually exclusive')
    if args.restart_from and args.seed_multiplier != 1.0:
        parser.error('--restart-from cannot be combined with --seed-multiplier')
    if ((args.restart_from or args.over_ionised_domain_guard or
         args.seed_multiplier != 1.0) and not args.development_unqualified):
        parser.error('seed/restart diagnostics are development-route only')

    run_directory = Path(args.run_directory).resolve()
    run_directory.mkdir(parents=True, exist_ok=True)
    if args.artifact:
        artifact = Path(args.artifact).resolve()
        manifest = json.loads((artifact / 'manifest.json').read_text())
        if file_hash(artifact / 'table.h5') != manifest['artifact_sha256']:
            parser.error('artifact table hash differs from manifest')
    else:
        spec_path = os.environ.get('EEDF_REAL_SPEC')
        if not spec_path:
            parser.error('EEDF_REAL_SPEC must name the pinned pure-Ar generation spec')
        artifact, manifest = generate_real_artifact(
            spec_path, run_directory / 'fixture', command_line=['runRealEEDFArgon.py'])
    restart_state = None
    if args.restart_from:
        restart_state = json.loads(Path(args.restart_from).read_text())
        if restart_state.get('scientific_status') != DEVELOPMENT_UNQUALIFIED_STATUS:
            parser.error('restart checkpoint lacks the failed-qualification development stamp')
        initialization = {
            'policy': 'restart from intermediate accepted development state',
            'provenance': str(Path(args.restart_from).resolve()),
            'artifact_sha256': manifest['artifact_sha256'],
            'checkpoint_time_s': restart_state['time_s'],
            'species_amounts_mol': restart_state['species_amounts_mol'],
            'u': restart_state['u'],
        }
    elif args.over_ionised_domain_guard:
        initialization = {
            'policy': 'OVER-IONISED INITIAL-STATE DOMAIN-GUARD TEST',
            'provenance': ('former pure-Ar reference seed retained only as a '
                           'named negative fixture'),
            'artifact_sha256': manifest['artifact_sha256'],
            'electron_mole_fraction': 1.e-6,
            'argon_ion_mole_fraction': 1.e-6,
        }
    else:
        initialization = power_consistent_initialization(
            artifact, manifest, multiplier=args.seed_multiplier)
    deck_path = run_directory / 'input.py'
    deck_path.write_text(pure_argon_deck(
        artifact, manifest, seed_multiplier=args.seed_multiplier,
        restart_state=restart_state,
        over_ionised=args.over_ionised_domain_guard))

    status = {
        'artifact': str(artifact),
        'artifact_sha256': manifest['artifact_sha256'],
        'accepted': manifest['accepted'],
        'branches': manifest['branches'],
        'branch_certification': manifest['branch_certification'],
        'blockers': reactor_readiness(manifest),
        'deck': str(deck_path),
        'initialization': initialization,
        'fingerprints': {
            'table': manifest['artifact_sha256'],
            'solver': git_fingerprint(Path(__file__).resolve().parents[3]),
            'deck': file_hash(deck_path),
            'configuration': configuration_fingerprint(args, deck_path),
            'table_solver': manifest.get('solver'),
        },
    }
    if args.development_unqualified:
        status.update(development_diagnostic(status['blockers']))
    status_path = run_directory / 'status.json'
    write_json(status_path, status)
    print(json.dumps(status, indent=2, sort_keys=True))
    if args.prepare_only:
        return 0
    if status['blockers'] and not args.development_unqualified:
        print('END_TO_END_BLOCKED: ' + '; '.join(status['blockers']))
        return 3

    rmg_directory = run_directory / 'rmg'
    rmg_directory.mkdir(parents=True, exist_ok=True)
    rmg = RMG(input_file=str(deck_path), output_directory=str(rmg_directory))
    if args.development_unqualified:
        rmg.initialize()
        reactor = rmg.reaction_systems[0]
        database_directory = (
            getattr(getattr(rmg, 'database', None), 'directory', None) or
            settings.get('database.directory'))
        if database_directory:
            status['fingerprints']['database'] = git_fingerprint(database_directory)
        else:
            status['fingerprints']['database'] = {
                'status': 'NOT EXECUTED',
                'reason': 'RMG database object did not expose its directory',
            }
        write_json(status_path, status)
        core_species = rmg.reaction_model.core.species
        core_reactions = development_argon_model(rmg.reaction_model.core.reactions)
        logging.info(
            '%s: starting pure-Ar production-solver entry with blockers: %s',
            DEVELOPMENT_UNQUALIFIED_STATUS, '; '.join(status['blockers']))
        reactor.configure_development_run(
            progress_interval_seconds=args.progress_interval_seconds,
            wall_budget_seconds=args.wall_budget_seconds)
        try:
            with development_unqualified_table_route():
                simulation = reactor.simulate(
                    core_species, core_reactions, [], [], [], [],
                    model_settings=rmg.model_settings_list[0],
                    simulator_settings=rmg.simulator_settings_list[0])
        except Exception as error:
            if args.over_ionised_domain_guard:
                negative = negative_fixture_result(reactor, error)
                write_json(run_directory / 'negative-fixture.json', negative)
                print(json.dumps(negative, indent=2, sort_keys=True))
                return 0 if negative['outcome'] == 'PASS' else 5
            failure = development_failure(reactor, status['blockers'], error)
            write_json(run_directory / 'failure.json', failure)
            logging.exception(
                '%s: pure-Ar run stopped before steady state',
                DEVELOPMENT_UNQUALIFIED_STATUS)
            return 4
        logging.info(
            '%s: solver returned at t=%g s; steady_state_reached=%r',
            DEVELOPMENT_UNQUALIFIED_STATUS, simulation[5],
            reactor.steady_state_reached)
    else:
        rmg.execute()
        reactor = rmg.reaction_systems[0]
    if args.over_ionised_domain_guard:
        negative = {
            'fixture': 'OVER-IONISED INITIAL-STATE DOMAIN-GUARD TEST',
            'outcome': 'FAIL',
            'reason': 'obsolete 1e-6 seed unexpectedly reached a terminal state',
        }
        negative.update(development_diagnostic(status['blockers']))
        write_json(run_directory / 'negative-fixture.json', negative)
        return 5
    endpoint_classification = candidate_endpoint_classification(reactor)
    budget = manifest_values(dict(reactor.energy_budget))
    coordinates = reactor._eedf_coordinates(reactor.y)
    result = {
        'arm': (args.arm_label or ('restart' if args.restart_from else
                                   '{0:g}x'.format(args.seed_multiplier))),
        'endpoint_classification': endpoint_classification,
        'converged_to_steady_state': bool(reactor.steady_state_reached),
        'final_time_s': float(reactor.t),
        'u': float(reactor.y[reactor.te_index]),
        'EN_Td': float(reactor.eedf_row.EN_Td),
        'mean_energy_eV': float(reactor.eedf_row.swarm['mean_energy_eV']),
        'Te_eff_K': float(reactor.Te.value_si),
        'electron_density_m^-3': float(budget['n_e']),
        'Ar4s_fraction': float(coordinates.get('Ar4s_total', 0.)),
        'energy_budget': budget,
        'A6a_relative': float(budget['A6a_relative']),
        'A6a_passed': bool(budget['A6a_passed']),
        'A6b_relative': float(budget['A6b_relative']),
        'A6b_passed': bool(budget['A6b_passed']),
    }
    if args.development_unqualified:
        trajectory = list(reactor.development_trajectory)
        if not trajectory:
            raise RuntimeError('development run returned without accepted-state trajectory')
        endpoint = trajectory[-1]
        result.update({
            'composition': endpoint['composition'],
            'wall_flux_mol_s': endpoint['wall_flux_mol_s'],
            'electron_power_partition_W': endpoint['electron_power_partition_W'],
            'convergence_time_s': float(reactor.t),
            'transient_extrema': dict(reactor.development_transient_extrema),
        })
        a6b = (dict(reactor.development_a6b_failure)
               if reactor.development_a6b_failure is not None else {
                   'check': 'A6b LoKI-B/table field-power consistency',
                   'outcome': 'PASS',
                   'value': float(budget['A6b_relative']),
                   'tolerance': float(budget['A6b_tolerance']),
                   'numerator': float(budget['A6b_numerator']),
                   'denominator': float(budget['A6b_denominator']),
                   'artifact_sha256': manifest['artifact_sha256'],
               })
        a_posteriori = a_posteriori_resolve(reactor, artifact, run_directory)
        run_manifest = reactor.eedf_run_manifest()
        run_manifest.update({
            'initialization': initialization,
            'fingerprints': status['fingerprints'],
            'A6b_outcome': a6b,
            'a_posteriori_LoKI_B': a_posteriori,
        })
        result.update({
            'A6b_outcome': a6b,
            'a_posteriori_LoKI_B': a_posteriori,
            'run_manifest': run_manifest,
        })
        result.update(development_diagnostic(status['blockers']))
        write_json(run_directory / 'trajectory.json', trajectory)
        write_json(run_directory / 'checkpoint.json',
                   select_restart_checkpoint(trajectory))
    else:
        result['run_manifest'] = reactor.eedf_run_manifest()
    write_json(run_directory / 'result.json', result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
