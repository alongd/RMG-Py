#!/home/alon/anaconda3/envs/rmg_env/bin/python
"""Reproduce the reference, execute the exact adoption sketch, and publish atomically.

Exit zero means the frozen report reproduced, including an honestly recorded
scientific rejection. --require-acceptance also requires the adoption test to pass;
--require-all-mutations also requires every substantive specified mutant to be killed.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import uuid
import xml.etree.ElementTree as ET

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
import checks
import mutation_tests

DEFAULT_OUTPUT=Path('/home/alon/runs/i030-met-kernel-reference/rework2/verifications')


def sha(data):
    return hashlib.sha256(data).hexdigest()


def load_oracle():
    spec=importlib.util.spec_from_file_location('verified_rouse_oracle',ROOT/'reference/run.py')
    oracle=importlib.util.module_from_spec(spec);spec.loader.exec_module(oracle)
    return oracle


def verify_historical_uncertainty_probe(oracle,cfg):
    """Reproduce the separately labelled A02 reanalysis from its actual commit."""
    probe=json.loads((ROOT/'uncertainty_probe.json').read_text())
    prefix=probe['source_commit']+':test/rmgpy/kmc/fixtures/i030_kernel_reference/'
    raw=subprocess.check_output(['git','show',prefix+'reference/results/results.json'],cwd=ROOT)
    candidate_raw=subprocess.check_output(['git','show',prefix+'candidate_results.json'],cwd=ROOT)
    assert sha(raw)==probe['A02_results_sha256']
    assert sha(candidate_raw)==probe['A02_candidate_report_sha256']
    old=json.loads(raw)
    old_rows={(row['case'],row['pair_class']):row for row in old['qualification']}
    new=oracle.qualify(old['simulation'],cfg,old['mapped'])
    candidates=json.loads(candidate_raw)['comparison']
    import math
    expected=[]
    for row,candidate in zip(new,candidates):
        previous=old_rows[row['case'],row['pair_class']]
        n=row['numerical_uncertainty']
        deviation=abs(math.log(candidate['candidate_over_reference']))
        expected.append({'case':row['case'],'pair_class':row['pair_class'],
            'A02_sum_marginal_log':previous['total_log_tolerance'],
            'A02_means_joint_method_log':n['numerical_log_band'],
            'observed_audit_contrasts_log':n['observed_audit_contrasts_log'],
            'joint_audit_sampling_log':n['joint_audit_sampling_log'],
            'primary_statistics_log':n['statistics_log'],
            'A02_candidate_over_reference':candidate['candidate_over_reference'],
            'A02_absolute_log_deviation':deviation,
            'A02_joint_numerical_consumed':min(deviation,n['numerical_log_band']),
            'A02_residual_model_log_required':max(0.,deviation-n['numerical_log_band']),
            'proposed_model_log':0.})
    assert expected==probe['rows'],'historical uncertainty probe changed'


def write_json(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def publish(stage, parent, identity):
    # Publication is one same-filesystem rename of a complete immutable directory.
    files={str(path.relative_to(stage)):sha(path.read_bytes()) for path in sorted(stage.rglob('*'))
           if path.is_file() and path.name!='COMPLETE.json'}
    write_json(stage/'COMPLETE.json',{'identity':identity,'files':files})
    destination=parent/('verified-'+identity)
    assert not destination.exists()
    stage.rename(destination)
    return destination


def assemble():
    text=(ROOT/'pack.md').read_text()
    sections=[('REPRODUCED', (ROOT/'reference/results/results.md').read_text()),
              ('CANDIDATE',checks.render_candidate(json.loads((ROOT/'candidate_results.json').read_text()))),
              ('MUTATION',mutation_tests.render_mutations(json.loads((ROOT/'mutation_results.json').read_text())))]
    for label,table in sections:
        pattern=r'<!-- BEGIN '+label+r' NUMBERS -->.*?<!-- END '+label+r' NUMBERS -->'
        text,count=re.subn(pattern,lambda match:table.rstrip(),text,flags=re.S)
        assert count==1, f'missing {label} marker'
    (ROOT/'pack.md').write_text(text)
    print('ASSEMBLED: all measured tables embedded in pack.md',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
    parser.add_argument('--workers',type=int,choices=range(1,9),default=8)
    parser.add_argument('--assemble',action='store_true')
    parser.add_argument('--mutations-only',action='store_true')
    parser.add_argument('--require-acceptance',action='store_true')
    parser.add_argument('--require-all-mutations',action='store_true')
    args=parser.parse_args()
    if args.assemble:
        assemble();return
    # Restrict this campaign and all descendants to the same eight allowed CPUs.
    if hasattr(os,'sched_getaffinity'):
        allowed=sorted(os.sched_getaffinity(0))
        os.sched_setaffinity(0,set(allowed[:8]))
    env=dict(os.environ,PYTEST_DISABLE_PLUGIN_AUTOLOAD='1',PYTHONDONTWRITEBYTECODE='1',
             MET_KERNEL_REFERENCE_WORKERS=str(args.workers))
    env.pop('MET_KERNEL_REFERENCE_PACK',None)  # Prove the ordinary fixture default.
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        env[key]='1'
    oracle=load_oracle()
    cfg,digest=oracle.read_parameters()
    program_paths=['reference/run.py','checks.py','audit_sampler.py','mutation_tests.py','verify_pack.py',
                   'proposed_adoption_test.py','guardTest.py','acceptance_policy.json','decision_order.json',
                   'uncertainty_probe.json']
    program_hashes={name:sha((ROOT/name).read_bytes()) for name in program_paths}
    pack=(ROOT/'pack.md').read_text()
    frozen=json.loads((ROOT/'reference/results/results.json').read_text())
    assert frozen['parameters_sha256']==digest,'frozen report uses different parsed inputs'
    versions={'python':oracle.platform.python_version(),'numpy':oracle.np.__version__,
              'scipy':oracle.scipy.__version__,'rng':'numpy.PCG64',
              'platform':oracle.platform.system(),'machine':oracle.platform.machine()}
    checks.assert_generation_history(cfg,frozen,versions)
    policy_commit=checks.assert_decision_order()
    policy_in_git=subprocess.check_output(['git','show',policy_commit+':test/rmgpy/kmc/fixtures/i030_kernel_reference/acceptance_policy.json'],cwd=ROOT)
    assert policy_in_git==(ROOT/'acceptance_policy.json').read_bytes(),'policy changed after its pre-measurement commit'
    subprocess.run(['git','merge-base','--is-ancestor',policy_commit,'HEAD'],cwd=ROOT,check=True)
    verify_historical_uncertainty_probe(oracle,cfg)
    args.output=args.output.resolve()
    args.output.mkdir(parents=True,exist_ok=True)
    identity=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')+'-'+uuid.uuid4().hex
    stage=Path(tempfile.mkdtemp(prefix='.pending-'+identity+'-',dir=args.output))
    if args.mutations_only:
        mutation=mutation_tests.run_mutations(cfg,frozen,artifact_dir=stage/'mutation-artifacts')
        write_json(stage/'verification.json',{'parameters_sha256':digest,'program_sha256':program_hashes,
                                               'versions':frozen['versions'],'pytest_version':__import__('pytest').__version__})
        write_json(stage/'mutation_results.json',mutation)
        (stage/'mutation_results.md').write_text(mutation_tests.render_mutations(mutation))
        destination=publish(stage,args.output,identity)
        print(f'MUTATION REPORT: {destination}',flush=True)
        print(mutation_tests.render_mutations(mutation),flush=True)
        if args.require_all_mutations and not mutation['required_target_mutations_all_killed']:
            print('MUTATION REQUIREMENT FAIL: a required substantive mutant survived.',flush=True)
            raise SystemExit(1)
        return
    match=re.search(r'```python\n(def test_met_kernel_reference\(tmp_path\):.*?)(?:\n```)',pack,re.S)
    assert match,'no executable adoption sketch in pack'
    test_path=stage/'test_met_kernel_reference.py'
    adoption_source=(ROOT/'proposed_adoption_test.py').read_text()
    assert adoption_source[adoption_source.index('def test_met_kernel_reference'):].strip()==match.group(1).strip(), 'pack differs from directly runnable proposed test'
    test_path.write_text(adoption_source)
    junit=stage/'adoption-junit.xml'
    subprocess.run([sys.executable,'-B','-m','pytest','-c','/dev/null','--rootdir',str(ROOT),
                    '--basetemp',str(stage/'guard-work'),'-p','no:cacheprovider',
                    '--junitxml',str(stage/'guard-junit.xml'),'-q',str(ROOT/'guardTest.py')],
                   cwd=ROOT,env=env,check=True)
    print('SKETCH VERIFIER: executing the exact proposed adoption test with fresh trajectories.',flush=True)
    completed=subprocess.run([sys.executable,'-B','-m','pytest','-c','/dev/null','--rootdir',str(ROOT),
                              '--basetemp',str(stage/'pytest-work'),'-p','no:cacheprovider',
                              '--junitxml',str(junit),'-q','-s',str(ROOT/'proposed_adoption_test.py')],cwd=ROOT,env=env)
    outputs=list((stage/'pytest-work').rglob('results.json'))
    assert len(outputs)==1,'expected one complete fresh oracle output'
    output=outputs[0]
    fresh=json.loads(output.read_text())
    checks.assert_generation_history(cfg,fresh,versions)
    assert {key:value for key,value in fresh.items() if key!='generation_provenance'}=={
        key:value for key,value in frozen.items() if key!='generation_provenance'
    },'fresh independent scientific report differs from the frozen report'
    assert fresh['generation_provenance']==[
        {'precision_stage':stage_index,'parameters_sha256':digest,
         'program_sha256':program_hashes['reference/run.py'],
         'driver_program_sha256':program_hashes['reference/run.py'],
         'cases':[case['name'] for case in cfg['cases']
                  if stage_index < case.get('precision_stages',cfg['precision_stages'])]}
        for stage_index in range(max(case.get('precision_stages',cfg['precision_stages']) for case in cfg['cases']))]
    rendered=oracle.markdown(fresh)
    assert rendered==output.with_name('results.md').read_text()
    assert rendered.rstrip() in pack,'fresh numeric table differs from pack.md'
    candidate=checks.target_checks(checks.load_target(),cfg,fresh)
    checks.assert_transport_and_branches(candidate)
    assert candidate==json.loads((ROOT/'candidate_results.json').read_text()),'candidate report changed'
    assert checks.render_candidate(candidate).rstrip() in pack
    # Distinguish a scientific rejection from a broken test or interrupted run.
    if candidate['adoption_pass']:
        assert completed.returncode==0,'accepted candidate failed adoption test'
    else:
        failures=ET.parse(junit).findall('.//testcase/failure')
        assert completed.returncode==1 and len(failures)==1
        text=(failures[0].text or '')+failures[0].get('message','')
        expected=('Rouse oracle is not numerically qualified' if not all(r['reference_qualified'] for r in candidate['comparison'])
                  else 'kernel exceeds numerical plus proposed model band')
        assert expected in text,'adoption sketch failed for an unexpected reason'
    mutation=mutation_tests.run_mutations(cfg,fresh,artifact_dir=stage/'mutation-artifacts')
    assert mutation==json.loads((ROOT/'mutation_results.json').read_text()),'mutation results changed'
    assert mutation_tests.render_mutations(mutation).rstrip() in pack
    spec=importlib.util.spec_from_file_location('live_sampler_audit',ROOT/'audit_sampler.py')
    audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
    audits=audit.run_audits(cfg)
    assert audits==json.loads((ROOT/'sampler_audits.json').read_text()),'independent sampler audits changed'
    adoption_audits=list((stage/'pytest-work').rglob('sampler_audits.json'))
    assert len(adoption_audits)==1 and json.loads(adoption_audits[0].read_text())==audits, 'adoption test did not execute every complete sampler audit'
    assert {name:sha((ROOT/name).read_bytes()) for name in program_paths}==program_hashes,'program changed during verification'
    # Parsed snapshot remains the fingerprint source; no post-simulation parameter rehash.
    write_json(stage/'verification.json',{'parameters_sha256':digest,'program_sha256':program_hashes,
                                         'met_source_sha256':candidate['met_source_sha256'],
                                         'versions':fresh['versions'],'pytest_version':__import__('pytest').__version__,
                                         'fresh_scientific_reference_equal':True,
                                         'authoring_generation_provenance':frozen['generation_provenance'],
                                         'fresh_generation_provenance':fresh['generation_provenance'],
                                         'pack_root_environment_unset':True,
                                         'complete_sampler_audits_inside_adoption_test':True,
                                         'model_tolerance_proposal':checks.read_policy(cfg)[0]['model_tolerance'],
                                         'policy_commit':policy_commit,
                                         'policy_commit_bytes_equal':True,
                                         'declared_before_corrected_measurement_and_comparison':True,
                                         'historical_uncertainty_probe_reproduced':True,
                                         'adoption_pass':candidate['adoption_pass'],
                                         'required_target_mutations_all_killed':mutation['required_target_mutations_all_killed']})
    for name,value in [('candidate_results.json',candidate),('mutation_results.json',mutation),('sampler_audits.json',audits)]:
        write_json(stage/name,value)
    shutil.copyfile(output,stage/'results.json')
    (stage/'results.md').write_text(rendered)
    shutil.copyfile(output.with_name('runtime.json'),stage/'runtime.json')
    destination=publish(stage,args.output,identity)
    print(f'REPRODUCTION PASS: all {len(fresh["simulation"])} reported ensembles, contact hashes, covariance samples, audits and tables reproduced exactly.',flush=True)
    print(f'PROPOSED SCIENTIFIC TEST: {"PASS" if candidate["adoption_pass"] else "FAIL"}; all required target mutations killed: {mutation["required_target_mutations_all_killed"]}',flush=True)
    print(f'ATOMIC PUBLICATION: {destination}',flush=True)
    if args.require_acceptance and not candidate['adoption_pass']:
        raise SystemExit(1)
    if args.require_all_mutations and not mutation['required_target_mutations_all_killed']:
        raise SystemExit(1)


if __name__=='__main__':
    main()
