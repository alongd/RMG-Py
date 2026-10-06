"""Summarize the frozen-table seed sweep and restart without changing results."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from real_eedf_harness import compare_development_endpoints
from rmgpy.solver.eedf_provider import DEVELOPMENT_UNQUALIFIED_STATUS


def build_summary(arms, restart):
    """Build the stamped sweep/restart comparison without changing results."""
    for result in list(arms) + [restart]:
        if result.get('converged_to_steady_state') is not True:
            raise ValueError(
                '{} did not converge to steady state'.format(
                    result.get('arm', '<unnamed arm>')))
    central = next(arm for arm in arms if arm['arm'] == '1x')
    summary = {
        'scientific_status': DEVELOPMENT_UNQUALIFIED_STATUS,
        'export_allowed': False,
        'qualification_allowed': False,
        'sweep': compare_development_endpoints(arms),
        'restart': compare_development_endpoints([central, restart]),
        'arms': arms,
        'restart_arm': restart,
    }
    summary['different_branches_detected'] = not (
        summary['sweep']['all_endpoints_agree'] and
        summary['restart']['all_endpoints_agree'])
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    parser.add_argument('--restart-result', required=True)
    parser.add_argument('arm_results', nargs='+')
    args = parser.parse_args()

    arms = [json.loads(Path(path).read_text()) for path in args.arm_results]
    labels = {arm['arm'] for arm in arms}
    required = {'0.1x', '1x', '10x', '20x'}
    if labels != required:
        parser.error('sweep arms must be exactly ' + repr(sorted(required)))
    restart = json.loads(Path(args.restart_result).read_text())
    summary = build_summary(arms, restart)
    Path(args.output).write_text(json.dumps(summary, indent=2, sort_keys=True))
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 2 if summary['different_branches_detected'] else 0


if __name__ == '__main__':
    raise SystemExit(main())
