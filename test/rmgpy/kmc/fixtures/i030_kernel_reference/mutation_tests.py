#!/home/alon/anaconda3/envs/rmg_env/bin/python
"""Run actual in-memory target mutations and independent sampler controls.

A surviving no-op is reported as a survivor, even if scientific adoption was
already rejected. Never count failure of an unmodified baseline as a kill.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import subprocess
import tempfile
import xml.etree.ElementTree as ET

for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
import checks


def replace_once(old,new):
    def transform(source):
        assert source.count(old)==1, f'mutation target not unique: {old!r}'
        return source.replace(old,new,1)
    return transform


def adoption_sampler_mutation(cfg, artifact_dir=None):
    """Execute the exact adoption test; N16-only excess OU noise fails preflight.

    A passing complete audit phase is the positive control. The unmodified full
    scientific comparison may fail, and is never used as a sampler-mutation kill.
    """
    spec=importlib.util.spec_from_file_location('adoption_audit_control',ROOT/'audit_sampler.py')
    audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
    baseline=audit.run_audits(cfg)
    source=(ROOT/'proposed_adoption_test.py').read_text()
    fixture='''import pytest
import sys
sys.path.insert(0, PACK_PATH)
@pytest.fixture(autouse=True)
def n16_ou_noise_mutation(monkeypatch):
    import audit_sampler
    import numpy as np
    original = audit_sampler.oracle.advance_modes
    coefficients = audit_sampler.oracle.mode_coefficients
    masks = {}
    def marked_coefficients(cfg, mapped, pair):
        rates, weights, variances = coefficients(cfg, mapped, pair)
        masks[id(rates)] = np.asarray([n == 16 for n in
            (mapped['units_i'], mapped['units_j']) for _ in range(n-1)])
        return rates, weights, variances
    def mutated(q, dt, rates, variances, rng):
        result = original(q, dt, rates, variances, rng)
        mask = masks[id(rates)]
        if mask.any():
            mean = q*np.exp(-dt[:,None]*rates)[:,:,None]
            result[:,mask,:] = mean[:,mask,:] + 2*(result[:,mask,:]-mean[:,mask,:])
        return result
    monkeypatch.setattr(audit_sampler.oracle, 'mode_coefficients', marked_coefficients)
    monkeypatch.setattr(audit_sampler.oracle, 'advance_modes', mutated)
'''.replace('PACK_PATH',repr(str(ROOT)))
    with tempfile.TemporaryDirectory(prefix='i030-adoption-mutant-') as temporary:
        directory=Path(temporary) if artifact_dir is None else Path(artifact_dir)
        directory.mkdir(parents=True,exist_ok=True)
        test=directory/'test_adoption_n16_noise.py'
        test.write_text(fixture+'\n'+source)
        junit=directory/'n16-noise-junit.xml'
        env=dict(os.environ,PYTEST_DISABLE_PLUGIN_AUTOLOAD='1',PYTHONDONTWRITEBYTECODE='1',
                 MET_KERNEL_REFERENCE_PACK=str(ROOT))
        result=subprocess.run([sys.executable,'-B','-m','pytest','-c','/dev/null',
            '--rootdir',str(directory),'-p','no:cacheprovider','--basetemp',str(directory/'pytest-work'),
            '--junitxml',str(junit),'-q',str(test)],env=env,capture_output=True,text=True)
        (directory/'pytest.log').write_text(result.stdout+result.stderr)
        failures=ET.parse(junit).findall('.//testcase/failure')
        messages=[(row.text or '')+row.get('message','') for row in failures]
        killed=(result.returncode==1 and len(failures)==1
                and 'actual sampled stationary variance failed' in messages[0]
                and not list((directory/'pytest-work').rglob('results.json')))
    return {'mutation':'double_OU_noise_N16_only_in_exact_adoption_test',
            'baseline_complete_audit_phase_pass':True,
            'baseline_audits_sha256':hashlib.sha256(json.dumps(baseline,sort_keys=True).encode()).hexdigest(),
            'adoption_test_source_sha256':hashlib.sha256(source.encode()).hexdigest(),
            'mutation_fixture_sha256':hashlib.sha256(fixture.encode()).hexdigest(),
            'pytest_returncode':result.returncode,'pytest_failure_count':len(failures),
            'killed':killed,'failure_stage':'stationary_variance' if killed else 'unexpected',
            'predicate':'only modes belonging to length16 chains, including either unequal orientation; length4/8 modes unchanged',
            'baseline_control_scope':'complete preflight; full unmodified scientific outcome reported separately'}


def run_mutations(cfg, reference, sampler=True, artifact_dir=None):
    baseline=checks.target_checks(checks.load_target(),cfg,reference)
    checks.assert_transport_and_branches(baseline)
    old='capture_radius = max(SIGMA_CONTACT, 2.0 * radius_gyration)'
    definitions=[
      ('capture_radius_x0.5',replace_once(old,'capture_radius = 0.5 * max(SIGMA_CONTACT, 2.0 * radius_gyration)')),
      ('capture_radius_x2',replace_once(old,'capture_radius = 2.0 * max(SIGMA_CONTACT, 2.0 * radius_gyration)')),
      ('candidate_diffusion_x4',replace_once('return 4.0 * math.pi * N_A * spin_factor * diffusivity * capture_radius',
                                            'return 4.0 * math.pi * N_A * spin_factor * (4.0 * diffusivity) * capture_radius')),
      ('shared_chain_diffusivity_x4',replace_once('return d0 / units\n','return 4.0 * d0 / units\n')),
      ('drop_one_chain_diffusion',replace_once('diffusivity = arm.chain_diffusivity(temperature, units_i) + arm.chain_diffusivity(\n        temperature, units_j\n    )',
                                             'diffusivity = arm.chain_diffusivity(temperature, units_i)')),
      ('length_min_to_max',replace_once('min(units_i, units_j) / 6.0','max(units_i, units_j) / 6.0')),
      ('length_min_to_i',replace_once('min(units_i, units_j) / 6.0','units_i / 6.0')),
      ('length_min_to_j',replace_once('min(units_i, units_j) / 6.0','units_j / 6.0')),
      ('length_min_to_mean',replace_once('min(units_i, units_j) / 6.0','((units_i + units_j) / 2.0) / 6.0')),
      ('drop_sigma0_floor',replace_once(old,'capture_radius = 2.0 * radius_gyration'))]
    prefactors={}
    for case in cfg['cases']:
        profile={row['pair_class']:row for row in reference['qualification'] if row['case']==case['name']}
        ratio=profile['end/end']['reference_reduced']/profile['mid/mid']['reference_reduced']
        prefactors[(case['units_i'],case['units_j'],'end/end')]=1/ratio
        prefactors[(case['units_i'],case['units_j'],'mid/mid')]=ratio
    line='return 4.0 * math.pi * N_A * spin_factor * diffusivity * capture_radius'
    definitions.append(('swap_reference_prefactors_in_candidate',replace_once(line,
        line+' * '+repr(prefactors)+'.get((units_i, units_j, pair_class), 1.0)')))
    n4_only={key:value for key,value in prefactors.items() if key[:2]==(4,4)}
    definitions.append(('swap_reference_prefactors_N4_only_diagnostic',replace_once(line,
        line+' * '+repr(n4_only)+'.get((units_i, units_j, pair_class), 1.0)')))
    noop_transform=replace_once('if pair_class not in PAIR_CLASSES:\n',
            'pair_class = {"end/end": "mid/mid", "mid/mid": "end/end"}.get(pair_class, pair_class)\n    if pair_class not in PAIR_CLASSES:\n')
    outcomes=[]
    for name, transform in definitions:
        report=checks.target_checks(checks.load_target(transform),cfg,reference)
        failed_anchors=[row['name'] for row in report['anchors'] if not row['pass']]
        failed_branches=[f"{row['case']}:{row['pair_class']}" for row in report['literal_kernel_branches'] if not row['pass']]
        killed=bool(failed_anchors or failed_branches)
        changed=report['anchors'] != baseline['anchors'] or report['comparison'] != baseline['comparison']
        new_rate_rejections=[f"{row['case']}:{row['pair_class']}"
                             for original,row in zip(baseline['comparison'],report['comparison'])
                             if original['reference_qualified'] and original['within_error_budget']
                             and not row['within_error_budget']]
        outcomes.append({'mutation':name,'baseline_control_pass':True,'killed':killed,
                         'changed_observable':changed,'failed_anchors':failed_anchors,
                         'failed_literal_branches':failed_branches,
                         'new_qualified_rate_rejections':new_rate_rejections,
                         'within_rouse_budget_count':sum(row['within_error_budget'] for row in report['comparison']),
                         'failed_qualified_rate_comparisons':[f"{row['case']}:{row['pair_class']}" for row in report['comparison'] if row['reference_qualified'] and not row['within_error_budget']]})
    noop_report=checks.target_checks(checks.load_target(noop_transform),cfg,reference)
    noop={'mutation':'swap_target_class_labels_only','changed_observable':noop_report['comparison']!=baseline['comparison'],
          'expected_survivor':noop_report['anchors']==baseline['anchors'] and noop_report['comparison']==baseline['comparison']}
    # An oracle-backed positive control tests tolerance discrimination separately
    # from the implementation's already-failing science and exact rule checks.
    budgets=reference['qualification']
    model_log=checks.read_policy(cfg)[0]['model_tolerance']['symmetric_log']
    control=[]
    for name,factor in [('radius_x0.5',.5),('radius_x2',2.),('diffusion_x4',4.),('drop_one_equal_chain',.5)]:
        equal_cases={case['name'] for case in cfg['cases'] if case['units_i']==case['units_j']}
        eligible=[row for row in budgets if name!='drop_one_equal_chain' or row['case'] in equal_cases]
        failures=[f"{row['case']}:{row['pair_class']}" for row in eligible
                  if row['qualified'] and abs(math.log(factor)) > row['total_log_tolerance']+model_log]
        control.append({'mutation':name,'positive_control_count':sum(row['qualified'] for row in eligible),
                        'killed_by_qualified_rate_budget':bool(failures),'failures':failures})
    # The reference prefactors are reversed in the candidate above, as the
    # review asked. Separately test an oracle-backed site's ordering witness.
    oracle_swaps=[]
    z=cfg['statistical_sigmas']
    for case in cfg['cases']:
        rows={row['pair_class']:row for row in budgets if row['case']==case['name']}
        end,mid=rows['end/end'],rows['mid/mid']
        separation=end['reference_reduced']-mid['reference_reduced']
        error=z*math.hypot(end['reference_se_reduced'],mid['reference_se_reduced'])
        oracle_swaps.append({'case':case['name'],'unmodified_ordering_control_pass':separation>error,
                             'end_minus_mid':separation,'statistical_ordering_bound':error,
                             'swapped_oracle_prefactors_rejected':separation>error})
    sampler_controls=[]
    if sampler:
        spec=importlib.util.spec_from_file_location('mutated_rouse_sampler',ROOT/'reference/run.py')
        oracle=importlib.util.module_from_spec(spec);spec.loader.exec_module(oracle)
        # Restrict the live checks to N16 to keep mutation qualification cheap.
        setup=dict(cfg['sampler'],dt_reduced=.01,lags=[0,1,2,5,10])
        own=dict(cfg,cases=[case for case in cfg['cases'] if case['name']=='N16'],pair_classes=['end/end'],sampler=setup)
        unmodified=oracle.sampler_checks(own)
        assert all(row['pass'] for row in unmodified)
        original=oracle.advance_modes
        import numpy as np
        def doubled_noise(q,dt,rates,variances,rng):
            mean=q*np.exp(-dt[:,None]*rates)[:,:,None]
            return mean+2*(original(q,dt,rates,variances,rng)-mean)
        audit_spec=importlib.util.spec_from_file_location('sampler_variance_mutation',ROOT/'audit_sampler.py')
        audit=importlib.util.module_from_spec(audit_spec);audit_spec.loader.exec_module(audit)
        audit.oracle=oracle
        audit.stationary_variance_checks(own)
        oracle.advance_modes=doubled_noise
        try:
            audit.stationary_variance_checks(own)
        except AssertionError as exc:
            killed=True;reason=str(exc)
        else:
            killed=False;reason='mutated sampler passed'
        sampler_controls.append({'mutation':'double_actual_OU_noise','baseline_control_pass':True,
                                 'killed':killed,'reason':reason})
        oracle.advance_modes=original
        projected=oracle.projected_site
        oracle.projected_site=lambda q,weights:2*projected(q,weights)
        try:
            oracle.sampler_checks(own)
        except AssertionError as exc:
            killed=True;reason=str(exc)
        else:
            killed=False;reason='mutated internal coordinate sampler passed'
        sampler_controls.append({'mutation':'double_sampled_internal_coordinates','baseline_control_pass':True,
                                 'killed':killed,'reason':reason})
        sampler_controls.append(adoption_sampler_mutation(cfg,artifact_dir))
    required_rates={'capture_radius_x0.5','capture_radius_x2','candidate_diffusion_x4',
                    'shared_chain_diffusivity_x4','drop_one_chain_diffusion',
                    'swap_reference_prefactors_in_candidate'}
    rate_mutations_rejected=all(row['new_qualified_rate_rejections'] for row in outcomes
                               if row['mutation'] in required_rates)
    return {'baseline_transport_and_rule_pass':True,
            'baseline_proposed_policy_comparison_pass':baseline['adoption_pass'],
            'baseline_qualified_rate_pass_count':sum(row['reference_qualified'] and row['within_error_budget'] for row in baseline['comparison']),
            'baseline_qualified_rate_failures':[f"{row['case']}:{row['pair_class']}" for row in baseline['comparison'] if row['reference_qualified'] and not row['within_error_budget']],
            'target_mutations':outcomes,
            'qualified_rate_positive_controls':control,'oracle_prefactor_swap_witnesses':oracle_swaps,
            'sampler_mutations':sampler_controls,'label_only_noop':noop,
            'required_rate_mutations_all_rejected':rate_mutations_rejected,
            'required_target_mutations_all_killed':all(row['killed'] for row in outcomes
                                                    if not row['mutation'].endswith('_diagnostic'))
                                                and rate_mutations_rejected
                                                and all(row['killed_by_qualified_rate_budget'] for row in control)
                                                and all(row['killed'] for row in sampler_controls),
            'limitation':'The required prefactor mutation rescales candidate end/end by C_mid/C_end and mid/mid by C_end/C_mid, using the independent reference. A bare target label swap is a separate observational no-op; it is reported as an expected survivor, never substituted for the requested prefactor perturbation.'}


def render_mutations(result):
    lines=['<!-- BEGIN MUTATION NUMBERS -->','',
           'Target baseline passes literal-input and kernel-rule controls before each mutation. Scientific adoption failure is never counted as a mutation kill.',
           f"Unmodified strict scientific comparison passed: {result['baseline_proposed_policy_comparison_pass']}; qualified baseline rate rows passing: {result['baseline_qualified_rate_pass_count']}/18.",'',
           '| Actual target mutation | Changes observable | Rejected | Failed anchors | Failed literal branches | Failed qualified rate checks | New rate rejections from passing baseline rows |',
           '|---|---|---|---:|---:|---:|---:|---:|']
    for row in result['target_mutations']:
        lines.append(f"| {row['mutation']} | {row['changed_observable']} | {row['killed']} | {len(row['failed_anchors'])} | {len(row['failed_literal_branches'])} | {len(row['failed_qualified_rate_comparisons'])} | {len(row['new_qualified_rate_rejections'])} |")
    lines+=['','| Rate control against qualified oracle | Positive controls | Rejected by rate budget |','|---|---:|---|']
    for row in result['qualified_rate_positive_controls']:
        lines.append(f"| {row['mutation']} | {row['positive_control_count']} | {row['killed_by_qualified_rate_budget']} |")
    lines+=['','| Oracle prefactor swap witness | Original end > mid by 4 SE | Swap rejected |','|---|---|---|']
    for row in result['oracle_prefactor_swap_witnesses']:
        lines.append(f"| {row['case']} | {row['unmodified_ordering_control_pass']} | {row['swapped_oracle_prefactors_rejected']} |")
    lines+=['','| Sampler mutation | Rejected | Failure stage |','|---|---|---|']
    for row in result['sampler_mutations']:
        lines.append(f"| {row['mutation']} | {row['killed']} | {row.get('failure_stage',row.get('reason',''))} |")
    lines+=['','The N16-only noise mutant executes the exact adoption-test body and must fail its stationary-variance audit before any trajectory generation. Its positive control is the complete unmodified audit phase, not an already-failing scientific comparison.','',result['limitation'],'','<!-- END MUTATION NUMBERS -->']
    return '\n'.join(lines)+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference',type=Path,default=ROOT/'reference/results/results.json')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--require-all',action='store_true',help='exit 1 if any specified substantive target mutant survives')
    args=parser.parse_args()
    spec=importlib.util.spec_from_file_location('rouse_inputs',ROOT/'reference/run.py')
    oracle=importlib.util.module_from_spec(spec);spec.loader.exec_module(oracle)
    cfg,digest=oracle.read_parameters()
    reference=json.loads(args.reference.read_text())
    assert reference['parameters_sha256']==digest
    artifact_dir=args.output.parent/(args.output.stem+'-artifacts') if args.output else None
    result=run_mutations(cfg,reference,artifact_dir=artifact_dir)
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(render_mutations(result),flush=True)
    if args.require_all and not result['required_target_mutations_all_killed']:
        print('MUTATION REQUIREMENT FAIL: a required substantive mutant survived.',flush=True)
        raise SystemExit(1)


if __name__=='__main__':
    main()
