"""Generate the real table, prepare the pure-Ar deck, and run it when admissible."""

import argparse
import json
import logging
import math
import os
from pathlib import Path
import sys
import traceback

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from rmgpy.rmg.main import RMG
from rmgpy.solver.eedf_provider import development_unqualified_table_route
from rmgpy.solver.electronegative import manifest_values
from rmgpy.tools.eedf.schema import file_hash
from real_eedf_harness import (
    generate_real_artifact,
    development_argon_model,
    pure_argon_deck,
    reactor_readiness,
)


def development_diagnostic(blockers):
    """Return the mandatory warning fields for every harness diagnostic."""
    return {
        'scientific_status': 'DEVELOPMENT \u2014 UNQUALIFIED TABLE',
        'artifact_blockers': list(blockers),
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
    args = parser.parse_args()
    if (not math.isfinite(args.progress_interval_seconds) or
            args.progress_interval_seconds <= 0.):
        parser.error('--progress-interval-seconds must be finite and positive')
    if args.wall_budget_seconds is not None:
        if not math.isfinite(args.wall_budget_seconds) or args.wall_budget_seconds <= 0.:
            parser.error('--wall-budget-seconds must be finite and positive')
        if not args.development_unqualified:
            parser.error('--wall-budget-seconds is development-route only')

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
    deck_path = run_directory / 'input.py'
    deck_path.write_text(pure_argon_deck(artifact, manifest))

    status = {
        'artifact': str(artifact),
        'artifact_sha256': manifest['artifact_sha256'],
        'accepted': manifest['accepted'],
        'branches': manifest['branches'],
        'branch_certification': manifest['branch_certification'],
        'blockers': reactor_readiness(manifest),
        'deck': str(deck_path),
    }
    if args.development_unqualified:
        status.update(development_diagnostic(status['blockers']))
    status_path = run_directory / 'status.json'
    status_path.write_text(json.dumps(status, indent=2, sort_keys=True))
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
        core_species = rmg.reaction_model.core.species
        core_reactions = development_argon_model(rmg.reaction_model.core.reactions)
        logging.info(
            'DEVELOPMENT \u2014 UNQUALIFIED TABLE: starting pure-Ar production-solver '
            'entry with blockers: %s', '; '.join(status['blockers']))
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
            failure = development_failure(reactor, status['blockers'], error)
            (run_directory / 'failure.json').write_text(
                json.dumps(failure, indent=2, sort_keys=True))
            logging.exception(
                'DEVELOPMENT \u2014 UNQUALIFIED TABLE: pure-Ar run stopped before '
                'steady state')
            return 4
        logging.info(
            'DEVELOPMENT \u2014 UNQUALIFIED TABLE: solver returned at t=%g s; '
            'steady_state_reached=%r', simulation[5], reactor.steady_state_reached)
    else:
        rmg.execute()
        reactor = rmg.reaction_systems[0]
    budget = manifest_values(dict(reactor.energy_budget))
    coordinates = reactor._eedf_coordinates(reactor.y)
    result = {
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
        'run_manifest': reactor.eedf_run_manifest(),
    }
    if args.development_unqualified:
        result.update(development_diagnostic(status['blockers']))
    (run_directory / 'result.json').write_text(json.dumps(result, indent=2, sort_keys=True))
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
