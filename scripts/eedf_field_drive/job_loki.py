#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""Job-only LoKI preparation and bounded failure evidence; solver files stay unchanged."""
import functools
import fcntl
import resource
import time
import json
import os
from pathlib import Path
import re
import subprocess
import shutil

from rmgpy.tools.eedf.integrity import _solver_environment
from rmgpy.tools.eedf.loki import LoKIDriver
from rmgpy.tools.eedf.schema import _file_hash

LOG_LIMIT = 65536


def prepare_spec(root):
    root = Path(root).resolve()
    results = root / 'results'
    results.mkdir(parents=True, exist_ok=True)
    spec = json.loads((root / 'spec-template.json').read_text())
    # Runtime scheduling only: keep frozen physical inputs and qualification criteria.
    # The solver uses a clean environment derived from spec, not inherited OMP settings.
    spec.update(omp_threads=int(os.environ['OMP_NUM_THREADS']))
    assert spec['omp_threads'] == 5
    assert spec['max_parallel'] == 3 and spec['timeout_s'] == 6200
    print('SOLVER_RESOURCES', spec['omp_threads'], spec['max_parallel'], flush=True)
    for item in spec['input_files'].values():
        item['path'] = str(root / item['path'])
        assert _file_hash(item['path']) == item['sha256']
    spec['channel_map']['path'] = str(root / spec['channel_map']['path'])
    for name, path in [('binary', root / 'loki-build/app/loki'),
                       ('cmake_cache', root / 'loki-build/CMakeCache.txt')]:
        spec[name] = {'path': str(path), 'sha256': _file_hash(path)}
    spec['compiler'] = subprocess.check_output([os.environ['CXX'], '--version'], text=True).splitlines()[0]
    linked = subprocess.run(['ldd', spec['binary']['path']], env=_solver_environment(spec),
                            capture_output=True, text=True, check=True)
    (results / 'loki-ldd.stdout.log').write_text(linked.stdout)
    (results / 'loki-ldd.stderr.log').write_text(linked.stderr)
    assert 'not found' not in linked.stdout and 'gcc14.1' not in linked.stdout, linked.stdout
    spec['shared_objects'] = {}
    for name in re.findall(r'(?:=>\s+|^\s*)(/\S+)', linked.stdout, re.M):
        path = Path(name).resolve()
        if path.is_relative_to(root / 'loki-build'):
            spec['shared_objects'][str(path)] = _file_hash(path)
    assert spec['shared_objects']
    spec.update(scratch_root=str(root / 'loki-runs'), output_root=str(root / 'tables'))
    record_build_provenance(root, spec)
    return spec


def capture_logs(root, work, category, reason):
    """Copy only the last 64 KiB per solver stream/setup, including partial files."""
    root, work = Path(root).resolve(), Path(work).resolve()
    if work.parent != root / 'loki-runs':
        raise ValueError('solver log path outside this job')
    target = root / 'results' / category / work.name
    target.mkdir(parents=True, exist_ok=True)
    record = {'job_name': work.name, 'reason': reason, 'limit_bytes': LOG_LIMIT,
              'files': {}, 'missing': []}
    for name in ('stdout.log', 'stderr.log', 'setup.in'):
        source = work / name
        if not source.is_file():
            record['missing'].append(name)
            continue
        with source.open('rb') as stream:
            stream.seek(0, 2)
            size = stream.tell()
            stream.seek(max(0, size - LOG_LIMIT))
            data = stream.read(LOG_LIMIT)
        (target / name).write_bytes(data)
        record['files'][name] = {'original_bytes': size, 'saved_bytes': len(data)}
    (target / 'capture.json').write_text(json.dumps(record, indent=2) + '\n')
    marker = ('LOKI_REFUSAL_CAPTURED' if category == 'loki-refusals'
              else 'LOKI_PREFLIGHT_LOGS')
    line = marker + ' ' + work.name + ' ' + str(target)
    print(line, flush=True)
    with (root / 'results/steps.log').open('a') as stream:
        stream.write(line + '\n')
    return target


def install_refusal_capture(root):
    """Preserve strict exceptions and return values; add evidence only."""
    original = LoKIDriver.run_prepared

    @functools.wraps(original)
    def recorded(self, bundle):
        try:
            return original(self, bundle)
        except Exception as error:
            try:
                work = Path(bundle.spec['scratch_root']) / bundle.job_name
                capture_logs(root, work, 'loki-refusals',
                             type(error).__name__ + ': ' + str(error))
            except Exception as capture_error:
                print('LOKI_LOG_CAPTURE_FAILED ' + str(capture_error), flush=True)
            raise
    LoKIDriver.run_prepared = recorded


def record_build_provenance(root, spec):
    """Bind the reviewed source patch to this build's actual normal solver pins."""
    approved = json.loads((root / 'growth-mixing-provenance.json').read_text())
    source_hash = _file_hash(root / 'loki-src/source/ElectronKinetics.cpp')
    if source_hash != approved['patched_source_sha256']:
        raise ValueError('LoKI build source differs from approved mixing patch')
    if _file_hash(root / 'growth-mixing.patch') != approved['patch_sha256']:
        raise ValueError('LoKI growth patch identity changed')
    record = {
        'scope': approved['scope'],
        'upstream_commit': approved['upstream_commit'],
        'patch_sha256': approved['patch_sha256'],
        'patched_source_sha256': source_hash,
        'compiler': spec['compiler'],
        'binary': spec['binary'],
        'cmake_cache': spec['cmake_cache'],
        'shared_objects': spec['shared_objects'],
        'validation_provenance_sha256': _file_hash(
            root / 'growth-mixing-provenance.json'),
    }
    for name in ('growth-mixing.patch', 'growth-mixing-provenance.json'):
        shutil.copyfile(root / name, root / 'loki-build' / name)
    text = json.dumps(record, indent=2) + '\n'
    (root / 'results/loki-build-provenance.json').write_text(text)
    (root / 'loki-build/growth-mixing-build.json').write_text(text)


def install_row_timing(root):
    """Record each isolated solver call, including refusals, in results/."""
    original = LoKIDriver.run_prepared
    target = Path(root) / 'results/row-timings.jsonl'
    target.parent.mkdir(parents=True, exist_ok=True)

    def emit(record):
        with target.open('a') as stream:
            fcntl.flock(stream, fcntl.LOCK_EX)
            stream.write(json.dumps(record, allow_nan=False) + '\n')
            stream.flush()
            fcntl.flock(stream, fcntl.LOCK_UN)

    @functools.wraps(original)
    def timed(self, bundle):
        record = {
            'job_name': bundle.job_name, 'pid': os.getpid(),
            'coordinate': bundle.coordinate_values,
            'fields': len(bundle.fields_Td), 'fields_Td': list(bundle.fields_Td),
            'timeout_s': bundle.spec['timeout_s'],
            'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        }
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        start = time.monotonic()
        emit(dict(record, event='start', status='running'))
        try:
            rows = original(self, bundle)
        except BaseException as error:
            record.update(status='refused', exception=type(error).__name__,
                          message=str(error))
            raise
        else:
            record.update(status='passed', returned_rows=len(rows))
            return rows
        finally:
            usage = resource.getrusage(resource.RUSAGE_CHILDREN)
            record.update(
                event='end', ended_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                wall_seconds=time.monotonic() - start,
                child_user_seconds=usage.ru_utime - before.ru_utime,
                child_system_seconds=usage.ru_stime - before.ru_stime)
            try:
                emit(record)
            except Exception:
                if record['status'] == 'passed':
                    raise
                # An evidence-write failure must not conceal the solver refusal.
                print('LOKI_TIMING_WRITE_FAILED ' + bundle.job_name, flush=True)
    LoKIDriver.run_prepared = timed
