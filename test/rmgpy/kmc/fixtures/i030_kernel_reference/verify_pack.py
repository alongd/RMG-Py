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

DEFAULT_OUTPUT=Path('/home/alon/runs/i030-met-kernel-reference/rework/verifications')


def sha(data):
    return hashlib.sha256(data).hexdigest()


def load_oracle():
    spec=importlib.util.spec_from_file_location('verified_rouse_oracle',ROOT/'reference/run.py')
    oracle=importlib.util.module_from_spec(spec);spec.loader.exec_module(oracle)
    return oracle


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
             MET_KERNEL_REFERENCE_PACK=str(ROOT),MET_KERNEL_REFERENCE_WORKERS=str(args.workers))
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        env[key]='1'
    oracle=load_oracle()
    cfg,digest=oracle.read_parameters()
    program_paths=['reference/run.py','checks.py','audit_sampler.py','mutation_tests.py','verify_pack.py',
                   'extend_floor.py','reference/provenance/measurement_04.py',
                   'reference/provenance/parameters_04.json']
    program_hashes={name:sha((ROOT/name).read_bytes()) for name in program_paths}
    pack=(ROOT/'pack.md').read_text()
    frozen=json.loads((ROOT/'reference/results/results.json').read_text())
    assert frozen['parameters_sha256']==digest,'frozen report uses different parsed inputs'
    versions={'python':oracle.platform.python_version(),'numpy':oracle.np.__version__,
              'scipy':oracle.scipy.__version__,'rng':'numpy.PCG64',
              'platform':oracle.platform.system(),'machine':oracle.platform.machine()}
    checks.assert_generation_history(cfg,frozen,versions)
    args.output=args.output.resolve()
    args.output.mkdir(parents=True,exist_ok=True)
    identity=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')+'-'+uuid.uuid4().hex
    stage=Path(tempfile.mkdtemp(prefix='.pending-'+identity+'-',dir=args.output))
    if args.mutations_only:
        mutation=mutation_tests.run_mutations(cfg,frozen)
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
    prelude='import json\nimport os\nfrom pathlib import Path\nimport subprocess\nimport sys\nsys.dont_write_bytecode=True\n'
    test_path.write_text(prelude+'\n'+match.group(1)+'\n')
    junit=stage/'adoption-junit.xml'
    print('SKETCH VERIFIER: executing the exact proposed adoption test with fresh trajectories.',flush=True)
    completed=subprocess.run([sys.executable,'-B','-m','pytest','-c','/dev/null','--rootdir',str(stage),
                              '--basetemp',str(stage/'pytest-work'),'-p','no:cacheprovider',
                              '--junitxml',str(junit),'-q','-s',str(test_path)],cwd=ROOT,env=env)
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
                  else 'kernel exceeds independent Rouse error budget')
        assert expected in text,'adoption sketch failed for an unexpected reason'
    mutation=mutation_tests.run_mutations(cfg,fresh)
    assert mutation==json.loads((ROOT/'mutation_results.json').read_text()),'mutation results changed'
    assert mutation_tests.render_mutations(mutation).rstrip() in pack
    spec=importlib.util.spec_from_file_location('live_sampler_audit',ROOT/'audit_sampler.py')
    audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
    audits=audit.run_audits(cfg)
    assert audits==json.loads((ROOT/'sampler_audits.json').read_text()),'independent sampler audits changed'
    assert {name:sha((ROOT/name).read_bytes()) for name in program_paths}==program_hashes,'program changed during verification'
    # Parsed snapshot remains the fingerprint source; no post-simulation parameter rehash.
    write_json(stage/'verification.json',{'parameters_sha256':digest,'program_sha256':program_hashes,
                                         'met_source_sha256':candidate['met_source_sha256'],
                                         'versions':fresh['versions'],'pytest_version':__import__('pytest').__version__,
                                         'fresh_scientific_reference_equal':True,
                                         'authoring_generation_provenance':frozen['generation_provenance'],
                                         'fresh_generation_provenance':fresh['generation_provenance'],
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
