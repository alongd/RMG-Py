#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the "Software"),  #
# to deal in the Software without restriction, including without limitation  #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""Reactor-level ownership and caching for fingerprinted EEDF tables."""

import copy
from collections.abc import Mapping
from contextlib import contextmanager
from contextvars import ContextVar
from functools import wraps
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import threading
import time
from types import MappingProxyType
import uuid
import weakref

import numpy as np

from rmgpy.solver.eedf import (
    AmbiguousBranch,
    EEDFError,
    EEDFRow,
    EEDFTable,
    FingerprintMismatch,
    IllConditionedCoordinate,
    OutOfDomain,
    _detached_copy,
)
from rmgpy.tools.eedf.schema import (
    _canonical_json,
    _content_hash,
    qualification_rule_for_quantity,
)
from rmgpy.tools.eedf.signal_guard import run_guarded_attempt
from rmgpy.solver._eedf_ownership import (
    Registry, ownership_scope, ownership_use, ownership_retire,
)


DEVELOPMENT_UNQUALIFIED_STATUS = (
    'DEVELOPMENT ONLY \u2014 TABLE QUALIFICATION FAILED')
DEVELOPMENT_DIAGNOSTIC_STATUS = 'DEVELOPMENT \u2014 not a qualification'
MANDATORY_QUALIFICATION_ARTIFACTS = (
    'setup.in',
    'interpolated-row.json',
    'solver-row.json',
    'qualification-evidence.json',
)
_development_unqualified_route = ContextVar(
    'eedf_development_unqualified_route', default=False)


class DevelopmentWallBudgetExceeded(RuntimeError):
    """The private development runner exhausted its configured wall budget."""


@contextmanager
def _development_unqualified_table_route():
    """Temporarily admit an unqualified table for an explicit development run.

    This process-local context is deliberately absent from the input DSL and
    reactor constructor state.  It therefore cannot be selected by ``rmg.py``,
    serialized into an input file, or restored through ``__reduce__``.
    """
    token = _development_unqualified_route.set(True)
    try:
        yield
    finally:
        _development_unqualified_route.reset(token)


def _unqualified_artifact_blockers(manifest):
    """Return the artifact's qualification blockers without changing them."""
    blockers = []
    verdicts = manifest.get('held_out_verdicts', [])
    failed = [verdict for verdict in verdicts if verdict.get('passed') is not True]
    if manifest.get('accepted') is not True:
        blockers.append('artifact manifest is not accepted')
    if failed:
        blockers.append(
            '{0} of {1} held-out verdicts failed'.format(len(failed), len(verdicts)))
    for branch, certification in sorted(manifest.get('branch_certification', {}).items()):
        if not str(certification).lower().startswith('certified'):
            blockers.append('{0}: {1}'.format(branch, certification))
    return tuple(blockers)


def _freeze(value):
    """Return a recursively immutable defensive value."""
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, np.ndarray):
        source = np.ascontiguousarray(value)
        # A read-only flag on an owning ndarray is reversible with
        # ``setflags(write=True)``.  Back the public row with immutable bytes so
        # callers cannot turn mutation back on and corrupt a cached evaluation.
        result = np.frombuffer(source.tobytes(), dtype=source.dtype).reshape(source.shape)
        return result
    return copy.deepcopy(value)


def _freeze_row(row):
    return EEDFRow(**{name: _freeze(value) for name, value in vars(row).items()})


def _qualification_locked(method):
    """Serialize one provider's qualification transaction and current-run latch."""
    @wraps(method)
    def locked(provider, *args, **kwargs):
        with provider._qualification_lock:
            return method(provider, *args, **kwargs)
    return locked


class _QualificationPublication:
    """One revocable publisher and the response held by its active frames."""

    def __init__(self):
        self.generation = {}
        self.record = {}
        self.commit = None
        self.record_path = None
        self.temporary = None
        self.file = None
        self.revoked = False
        self.published = False


class _QualificationCapture:
    """One registered private capture; only this context writes its namespace."""

    def __init__(self, owner):
        self.owner = owner
        self._resources = owner._registry.track(self, owner._resources)
        self.name = '.eedf-{0}.private'.format(uuid.uuid4().hex)
        self.directory = owner.private_parent / self.name

    @property
    def fd(self):
        return self._resources.descriptor(0)

    @ownership_use
    def capture(self):
        parent = self.owner.parent_fd
        try:
            os.mkdir(self.name, 0o700, dir_fd=self.owner.private_parent_fd)
            self._resources.open(0, self.name, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                                 dir_fd=self.owner.private_parent_fd)
            try:
                os.rename(self.owner.destination.name, 'captured',
                          src_dir_fd=parent, dst_dir_fd=self.fd)
            except FileNotFoundError:
                pass
        finally:
            self.settle()

    @ownership_use
    def settle(self):
        """Reconcile the actual captured inode, including interrupted syscalls."""
        try:
            if self.fd is None:
                return
            try:
                info = os.stat('captured', dir_fd=self.fd, follow_symlinks=False)
            except FileNotFoundError:
                return
            if (info.st_dev, info.st_ino) == self.owner.identity:
                os.unlink('captured', dir_fd=self.fd)
            else:
                self.owner.ambiguous = True
                try:
                    os.link('captured', self.owner.destination.name,
                            src_dir_fd=self.fd, dst_dir_fd=self.owner.parent_fd,
                            follow_symlinks=False)
                except OSError as error:
                    raise EEDFError(
                        'foreign file retained at {0}; original path {1}: {2}'.format(
                            self.directory / 'captured', self.owner.destination, error)) from error
                os.unlink('captured', dir_fd=self.fd)
        finally:
            # Never recursively remove a directory: it may hold foreign data.
            try:
                os.rmdir(self.name, dir_fd=self.owner.private_parent_fd)
            except FileNotFoundError:
                pass
            except OSError as error:
                raise EEDFError('private capture retained at ' + str(self.directory)) from error

class _QualificationExport:
    """Descriptor and private staging owned by a suspended single-file commit."""

    def __init__(self, destination, run_id=None, kind='export', *, registry):
        self._registry = registry
        self._resources = registry.track(self)
        self.destination = destination
        self.private_parent = destination.parent
        self.directory = destination.parent / ('.eedf-{0}.private'.format(uuid.uuid4().hex))
        self.source = self.directory / ('record.' + kind + '.tmp')
        self.identity = None
        self.captures = self._resources.captures
        self.ambiguous = False
        self.revoked = False
        self.run_id = run_id
        self.failure_payload = None

    @property
    def parent_fd(self):
        return self._resources.descriptor(0)

    @property
    def private_parent_fd(self):
        return self._resources.descriptor(1)

    @property
    def directory_fd(self):
        return self._resources.descriptor(2)

    @property
    def stream(self):
        return None if self._resources.closed else self._resources

    @ownership_use
    def anchor(self):
        """Keep staging outside the public directory when its parent is on this device."""
        self._resources.open(0, self.destination.parent,
                             os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        parent = self._resources.open(1, '..', os.O_RDONLY | os.O_DIRECTORY,
                                      dir_fd=self.parent_fd)
        if os.fstat(parent).st_dev == os.fstat(self.parent_fd).st_dev:
            self.private_parent = self.destination.parent.parent
        else:
            self._resources.duplicate(1, 0)
        self.directory = self.private_parent / self.directory.name
        self.source = self.directory / self.source.name

    @ownership_use
    def stage(self, payload):
        if self.parent_fd is None:
            self.anchor()
        os.mkdir(self.directory.name, 0o700, dir_fd=self.private_parent_fd)
        self._resources.open(2, self.directory.name,
                             os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                             dir_fd=self.private_parent_fd)
        fd = self._resources.open(3, self.source.name,
                                  os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                                  0o600, dir_fd=self.directory_fd)
        info = os.fstat(fd)
        self.identity = (info.st_dev, info.st_ino)
        if self.revoked:
            raise EEDFError('revoked export refuses staging')
        if self.stream.write(payload) != len(payload):
            raise EEDFError('incomplete export staging write')

    @ownership_use
    def _remove_owned(self, path):
        if self.identity is None or self.parent_fd is None or self.ambiguous:
            return
        try:
            os.stat(self.destination.name, dir_fd=self.parent_fd, follow_symlinks=False)
        except FileNotFoundError:
            return
        capture = _QualificationCapture(self)
        capture.capture()
        if self.ambiguous:
            raise EEDFError('export path ownership is ambiguous: ' + str(path))

    @ownership_use
    def finish(self):
        """Retire only the private alias, leaving a successful public inode intact."""
        if self.directory_fd is not None:
            try:
                info = os.stat(self.source.name, dir_fd=self.directory_fd, follow_symlinks=False)
            except FileNotFoundError:
                pass
            else:
                if self.identity is not None and (info.st_dev, info.st_ino) != self.identity:
                    raise EEDFError('foreign staging file retained at ' + str(self.source))
                os.unlink(self.source.name, dir_fd=self.directory_fd)
        if self.private_parent_fd is not None:
            try:
                os.rmdir(self.directory.name, dir_fd=self.private_parent_fd)
            except FileNotFoundError:
                pass

    def _settle_revocation(self):
        errors = []
        invalidation_error = None
        if self.identity is not None:
            try:
                marker = (self.failure_payload or
                          json.dumps({'status': 'FAIL', 'run_id': self.run_id}).encode())
                os.ftruncate(self.stream.fileno(), 0)
                if os.pwrite(self.stream.fileno(), marker, 0) != len(marker):
                    raise EEDFError('incomplete revocation marker')
            except BaseException as error:
                invalidation_error = error
        for capture in tuple(self.captures):
            try:
                capture.settle()
            except BaseException as error:
                errors.append(error)
        try:
            self._remove_owned(self.destination)
        except BaseException as error:
            errors.append(error)
        try:
            self.finish()
        except BaseException as error:
            errors.append(error)
        if errors:
            prefix = ('unlink and inode invalidation both failed: ' if invalidation_error
                      else 'inode invalidated but unlink failed: ')
            raise EEDFError(prefix + '; '.join(str(e) for e in errors)) from errors[-1]

    @ownership_retire
    def revoke(self, payload=None):
        self.revoked = True
        if payload is not None:
            self.failure_payload = payload
        completed = False
        try:
            self._settle_revocation()
            completed = True
        finally:
            if not completed:
                self._settle_revocation()

class _QualificationAdmission:
    """Revocable evidence; only commit_bytes consumes export authority.

    commit_bytes(destination, payload) accepts one new destination and bytes,
    never a callback. It stages a sibling file and atomically links its inode
    into the absent destination; existing files are never overwritten. Failure
    or a newer run revokes the capability and any suspended commit. Filesystems
    must support hard links. A path, record, manifest or boolean is not authority.
    """

    def __init__(self, provider, generation, commit, operation):
        self._provider = provider
        self._generation = generation
        self._commit = commit
        self._operation = operation

    def __bool__(self):
        raise TypeError('qualification admission is consumed only by commit_bytes')

    def _check(self):
        provider = self._provider
        if (provider._qualification_failure
                or provider._qualification_generation is not self._generation
                or not self._generation.get('published')
                or provider._qualification_commit is not self._commit):
            raise EEDFError('revoked qualification admission refuses ' + self._operation)

    @ownership_scope
    def commit_bytes(self, destination, payload):
        """Commit caller-prepared bytes to one absent path, or refuse.

        Revocation removes only this commit's exact path/inode. Ambiguous
        ownership refuses without deleting the replacement. If both unlink
        and descriptor-based invalidation fail, refusal reports the surviving
        file explicitly; filesystem permissions cannot be overridden.
        """
        if not isinstance(payload, bytes):
            raise TypeError('commit_bytes requires bytes, not a callback or mutable buffer')
        destination = Path(os.path.abspath(os.fspath(destination)))
        provider = self._provider
        with provider._qualification_lock:
            pending = _QualificationExport(
                destination, self._commit['run_id'], registry=provider._file_registry)
            provider._file_registry.register_export(pending, pending._resources)
            try:
                self._check()
                # Validate the exact immutable evidence again at consumption.
                provider._validate_acceptance_manifest(
                    self._commit['record_path'], self._commit['run_id'], self._commit)
                self._check()
                pending.stage(payload)
                self._check()
                info = pending.source.lstat()
                if (info.st_dev, info.st_ino) != pending.identity:
                    raise EEDFError('export source ownership is ambiguous')
                # Exclusive link publication cannot clobber another inode.
                # Reentry after this check removes/invalidates this very source.
                os.link(pending.source, destination, follow_symlinks=False)
                self._check()
                pending.finish()
                self._check()
                return destination
            except BaseException:
                pending.revoke()
                raise
            # Native scope exit follows this body's actual return, including
            # the return trace; then it closes and deregisters explicitly.


def _qualification_attempt_locked(method):
    """Hold the attempt lock and defer owned signals through FAILED finalization."""
    @wraps(method)
    def locked(provider, y, t, *args, **kwargs):
        context = None

        def attempt():
            nonlocal context
            context = provider._begin_qualification_publication()
            return method(provider, y, t, *args, **kwargs)

        def finalize(error):
            record = (provider._qualification_record
                      if provider._qualification_status == 'FAILED' else None)
            if record is not None:
                error = record.get('error', error)
            provider._enter_failed(error, y=y, record=record)
            return provider._qualification_record

        def complete(record):
            if context is None or not provider._publication_is_current(context):
                return provider._settled_qualification_record(y)
            try:
                provider._publish_acceptance(
                    context.commit, context.record_path)
            except BaseException as error:
                if not provider._publication_is_current(context):
                    if not isinstance(error, Exception):
                        raise
                    return provider._settled_qualification_record(y)
                publication_error = provider._enter_failed(
                    error, y=y, record=record)
                if not isinstance(error, Exception):
                    raise
                if publication_error is not None:
                    if not isinstance(publication_error, Exception):
                        raise publication_error
                    error = publication_error
                from rmgpy.exceptions import PlasmaStateError
                raise PlasmaStateError(
                    'EEDF qualification failed: ' + str(error)) from error
            return provider._settled_qualification_record(y)

        return run_guarded_attempt(
            provider._qualification_lock, attempt, finalize, complete)
    return locked


class EEDFProvider:
    """Immutable reactor-owned view of one accepted table branch.

    Rows returned by this provider are immutable handles issued by this exact
    provider. Consumers cannot pass a row from another table or provider into
    rate, transport, or energy accessors.
    """

    def __init__(self, path, current_model, *, artifact_sha256, reaction_map,
                 branch, empirical_laws=None):
        (path, current_model, artifact_sha256, reaction_map, branch,
         empirical_laws) = copy.deepcopy(
            (path, current_model, artifact_sha256, reaction_map, branch,
             empirical_laws))
        if not isinstance(artifact_sha256, str):
            raise TypeError('artifact_sha256 must be a string')
        if not isinstance(branch, str):
            raise TypeError('branch must be a string')
        development = bool(_development_unqualified_route.get())
        table = EEDFTable.load(path, current_model,
                               artifact_sha256=artifact_sha256,
                               require_accepted=not development)
        if not development and not table.a6b_source['passed']:
            raise FingerprintMismatch('A6b-source native field-power identity')
        if not development and not table.t10_interpolation['passed']:
            raise FingerprintMismatch('T10 interpolation qualification')
        if branch not in table.manifest['branches']:
            raise AmbiguousBranch('declared branch ' + str(branch))
        energy_index = table._layout[branch]['swarm/mean_energy_eV'][0]
        if np.any(np.diff(table._data[branch][..., energy_index], axis=0) <= 0.0):
            raise IllConditionedCoordinate(
                'mean energy must increase strictly with u=ln(E/N)')
        self._table = table
        # Reactor callbacks need repeated read-only access to artifact metadata.
        # Freeze one private snapshot here so those callbacks never pay for the
        # public manifest property's deliberately defensive deep copy.
        self._runtime_metadata = _freeze(table.manifest)
        self._branch = str(branch)
        self._reaction_map = MappingProxyType(self._validate_reaction_map(reaction_map))
        if empirical_laws is not None and not isinstance(empirical_laws, dict):
            raise TypeError('empirical_laws must be a mapping')
        self._empirical_laws = _freeze(empirical_laws or {})
        if development and table.manifest.get('accepted') is True:
            raise EEDFError(
                'development unqualified-table route refuses an already qualified artifact')
        self._development_blockers = (_unqualified_artifact_blockers(table.manifest)
                                      if development else ())
        self._cache_identity = None
        self._cached_row = None
        self._issued_rows = weakref.WeakValueDictionary()
        self._cache_lock = threading.RLock()
        self._qualification_lock = threading.RLock()
        self._qualification_phase = None
        self._qualification_failure = {}
        self._qualification_current = None
        self._qualification_generation = None
        self._qualification_publications = weakref.WeakSet()
        self._file_registry = Registry()
        self._qualification_exports = self._file_registry.pending_exports
        self._qualification_files = self._file_registry.files
        self._qualification_record = None
        self._qualification_record_path = None
        self._qualification_state_resolver = None
        self._qualification_run_id = None
        self._acceptance_manifest_path = None
        self._qualification_commit = None
        self._qualification_attempt = None

    @property
    def _qualification_status(self):
        """Failure dominates even a success marker written by a suspended call."""
        if self._qualification_failure:
            return 'FAILED'
        context = self._qualification_current
        if context is not None and context.published and not context.revoked:
            return 'PASS'
        return self._qualification_phase

    def _begin_qualification_publication(self):
        context = _QualificationPublication()
        self._qualification_generation = context.generation
        self._qualification_publications.add(context)
        self._qualification_current = context
        self._qualification_phase = None
        errors = self._revoke_publications(except_context=context)
        if errors:
            raise errors[0]
        return context

    def _publication_is_current(self, context):
        return (context is self._qualification_current
                and not context.revoked and not self._qualification_failure)

    def _revoke_publications(self, except_context=None):
        """Withdraw the actual rename sources of suspended publishers."""
        errors = []
        for pending in tuple(self._qualification_exports):
            try:
                pending.revoke()
            except BaseException as error:
                errors.append(error)
        for context in tuple(self._qualification_publications):
            if context is except_context:
                continue
            context.revoked = True
            if context.file is not None:
                try:
                    failure = self._qualification_failure.get('record', context.record)
                    context.file.revoke(
                        self._json_bytes({'status': 'FAIL',
                                          'run_id': failure.get('run_id')}))
                except BaseException as error:
                    errors.append(error)
        return errors

    def _sync_qualification_responses(self, record):
        """Update responses still held by outer completion/guard frames.

        Weak membership follows the actual frame lifetime, including its return;
        no explicit deregistration creates a gap before RETURN_VALUE. Completed
        independent runs retain their own snapshots once their contexts expire.
        """
        record = self._qualification_failure.get('record', record)
        for context in tuple(self._qualification_publications):
            context.record.clear()
            context.record.update(record)
        current = self._qualification_current
        self._qualification_record = current.record if current is not None else record
        return self._qualification_record

    def _settled_qualification_record(self, y):
        """Finish a failed publisher after a returning handler releases its I/O fault."""
        record = self._qualification_failure.get('record')
        if record is not None and self._qualification_record_path is None:
            self._enter_failed(record.get('error'), y=y, record=record)
        return self._qualification_record

    def _validate_reaction_map(self, reaction_map):
        if not isinstance(reaction_map, dict):
            raise TypeError('reaction_map must be a mapping')
        channels = self._runtime_metadata['channel_map']
        validated = {}
        for reaction_index, location in reaction_map.items():
            if (isinstance(reaction_index, bool) or
                    not isinstance(reaction_index, (int, np.integer)) or
                    reaction_index < 0):
                raise ValueError('invalid reaction index ' + repr(reaction_index))
            if not isinstance(location, (tuple, list)) or len(location) != 2:
                raise ValueError('reaction map entry must be (column, side)')
            column, side = location
            if (isinstance(column, bool) or
                    not isinstance(column, (int, np.integer)) or
                    not 0 <= column < len(channels)):
                raise ValueError('invalid channel column ' + repr(column))
            channel = channels[int(column)]
            if channel.get('classification') != 'A' or channel.get('reaction') is None:
                raise ValueError('channel column is not reaction-owned A chemistry')
            if side not in ('ine', 'sup'):
                raise ValueError('invalid channel side ' + repr(side))
            validated[int(reaction_index)] = (int(column), side)
        return validated

    @property
    def _branch_value(self):
        return self._branch

    @property
    def _reaction_mapping(self):
        return self._reaction_map

    @property
    def _empirical_law_metadata(self):
        return self._empirical_laws

    @property
    def _manifest_snapshot(self):
        """A defensive artifact metadata snapshot, without export authority."""
        manifest = copy.deepcopy(self._table.manifest)
        notice = self.development_notice
        if notice is not None:
            manifest.update(notice)
        return manifest

    @property
    def runtime_metadata(self):
        """Immutable callback metadata, without terminal export authority."""
        return self._runtime_metadata

    @property
    def development_notice(self):
        """Serializable warning carried by every development diagnostic."""
        if not self._development_blockers:
            return None
        return {
            'scientific_status': DEVELOPMENT_UNQUALIFIED_STATUS,
            'artifact_blockers': list(self._development_blockers),
            'export_allowed': False,
            'qualification_allowed': False,
        }

    @_qualification_locked
    def require_qualified(self, operation):
        """Return revocable current-run evidence for commit_bytes(path, bytes).

        Validation-time failure raises. Later failure may revoke an object
        already on the return stack; only its commit_bytes method can authorize
        terminal-qualified output. Records, statuses, paths and manifests are
        non-authoritative snapshots. Existing destinations are refused.
        """
        if not isinstance(operation, str):
            raise TypeError('operation must be a string')
        operation = str(operation)
        if self._development_blockers:
            raise EEDFError(
                '{0} refuses {1}: {2}'.format(
                    DEVELOPMENT_UNQUALIFIED_STATUS, operation,
                    '; '.join(self._development_blockers)))
        generation = self._qualification_generation
        if (self._qualification_failure or generation is None
                or not generation.get('published') or self._qualification_commit is None):
            status = self._qualification_status or 'NOT QUALIFIED'
            raise EEDFError('EEDF qualification {0} refuses {1}'.format(status, operation))
        admission = _QualificationAdmission(
            self, generation, self._qualification_commit, operation)
        try:
            self._validate_acceptance_manifest(
                self._acceptance_manifest_path, self._qualification_run_id,
                admission._commit)
            admission._check()
        except Exception as error:
            if (not self._qualification_failure
                    and self._qualification_generation is admission._generation):
                self._enter_failed(
                    error,
                    record=(dict(self._qualification_record)
                            if isinstance(self._qualification_record, Mapping)
                            else None),
                    record_path=self._acceptance_manifest_path,
                    run_id=self._qualification_run_id)
            raise EEDFError(
                'EEDF qualification manifest refuses {0}: {1}'.format(
                    operation, error)) from error
        return admission

    @property
    def qualification_record_path(self):
        """Diagnostic record locator; the path itself grants no export authority."""
        return self._qualification_record_path

    def bind_qualification_state(self, resolver):
        """Bind the reactor callback that describes an accepted terminal state."""
        if not callable(resolver):
            raise TypeError('qualification state resolver must be callable')
        self._qualification_state_resolver = resolver

    @staticmethod
    def _json_value(value):
        if isinstance(value, MappingProxyType):
            value = dict(value)
        if isinstance(value, dict):
            return {str(key): EEDFProvider._json_value(item)
                    for key, item in value.items()}
        if isinstance(value, (tuple, list)):
            return [EEDFProvider._json_value(item) for item in value]
        if isinstance(value, np.ndarray):
            return EEDFProvider._json_value(value.tolist())
        if isinstance(value, np.generic):
            return EEDFProvider._json_value(value.item())
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, float) and not np.isfinite(value):
            if np.isnan(value):
                return 'NaN'
            return 'Infinity' if value > 0. else '-Infinity'
        return value

    def _qualification_state(self, y, t):
        if isinstance(y, dict):
            state = copy.deepcopy(y)
        elif self._qualification_state_resolver is not None:
            state = copy.deepcopy(self._qualification_state_resolver(y, t))
        else:
            raise EEDFError('terminal qualification state is unavailable')
        self._capture_qualification_state(state)
        for name in (
                'u', 'composition', 'gas_fractions', 'state_populations',
                'state_statistical_weights', 'Tg_K', 'P_Pa'):
            if name not in state:
                if name in ('gas_fractions', 'state_populations',
                            'state_statistical_weights'):
                    raise EEDFError(
                        'terminal ' + name.replace('_', ' ') + ' are incomplete')
                raise EEDFError('terminal qualification state has no ' + name)
        state['t'] = float(t)
        state['u'] = float(state['u'])
        state['composition'] = {str(key): float(value)
                                for key, value in state['composition'].items()}
        if (not np.isfinite(state['t']) or not np.isfinite(state['u']) or
                not all(np.isfinite(value) for value in state['composition'].values())):
            raise EEDFError('terminal qualification state is non-finite')
        for name in ('Tg_K', 'P_Pa'):
            state[name] = float(state[name])
            if not np.isfinite(state[name]) or state[name] <= 0.:
                raise EEDFError('terminal qualification state has invalid ' + name)
        if 'scratch_root' in state:
            if not isinstance(state['scratch_root'], (str, os.PathLike)):
                raise EEDFError('terminal qualification state has invalid scratch_root')
            state['scratch_root'] = str(state['scratch_root'])
        return state

    def _capture_qualification_state(self, state):
        """Register a resolved output path before validating the state."""
        attempt = self._qualification_attempt
        if attempt is None:
            return
        if 'record_path' in state and state['record_path'] is not None:
            record_path = Path(state['record_path'])
            invalid = {
                'status': 'FAIL',
                'run_id': attempt['record'].get('run_id'),
                'comparison_passed': False,
                'error': 'qualification attempt has not committed PASS',
            }
            self._invalidate_acceptance(record_path, invalid)
            attempt['record_path'] = record_path
        attempt['record']['terminal_state'] = state

    @staticmethod
    def _write_json(path, value):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(EEDFProvider._json_value(value), indent=2,
                                   sort_keys=True, allow_nan=False) + '\n')
        return path

    @staticmethod
    def _json_bytes(value):
        return (json.dumps(EEDFProvider._json_value(value), indent=2,
                           sort_keys=True, allow_nan=False) + '\n').encode('utf-8')

    def _record_path(self, state, run_id):
        requested = state.get('record_path')
        if requested is not None:
            return Path(requested)
        scratch = state.get('scratch_root')
        if scratch is None:
            raise EEDFError('terminal qualification state has no scratch_root')
        return Path(scratch) / ('qualification-' + run_id) / 'acceptance.json'

    def _write_failure(self, record, record_path, run_id):
        record['status'] = 'FAIL'
        failure_path = record_path.with_name(
            '.{0}.{1}.failure.json'.format(record_path.stem, run_id))
        self._qualification_record = record
        written_record = record
        try:
            payload = self._json_bytes(written_record)
        except (TypeError, ValueError):
            written_record = {'status': 'FAIL', 'run_id': run_id,
                              'comparison_passed': False,
                              'error': 'qualification failure record was not serializable'}
            payload = self._json_bytes(written_record)
        for target in (failure_path, record_path):
            self._publish_failure(target, run_id, written_record, payload)
        self._qualification_record_path = record_path
        return self._qualification_record_path

    @ownership_scope
    def _publish_failure(self, target, run_id, record, payload):
        self._invalidate_record(target, payload)
        owned = self._new_record_file(target, run_id, 'failure')
        try:
            self._stage_record(owned, record)
            os.link(owned.source, owned.destination, follow_symlinks=False)
            owned.finish()
        except BaseException:
            owned.revoke(payload)
            raise

    def _new_record_file(self, path, run_id, kind):
        path = Path(os.path.abspath(os.fspath(path)))
        path.parent.mkdir(parents=True, exist_ok=True)
        owned = _QualificationExport(path, run_id, kind, registry=self._file_registry)
        self._file_registry.register_record(str(path), owned, owned._resources)
        owned.anchor()
        return owned

    def _stage_record(self, owned, record):
        """Stage bytes only in the registered process-private namespace."""
        owned.stage(self._json_bytes(record))

    @ownership_scope
    def _invalidate_record(self, path, payload):
        """Retire the pinned record inode; never reopen a replacement to write it."""
        path = Path(os.path.abspath(os.fspath(path)))
        owned = self._qualification_files.get(str(path))
        if owned is None:
            # An explicitly designated prior-process record is adopted only
            # for retirement, never as evidence that can issue admission.
            try:
                owned = _QualificationExport(path, registry=self._file_registry)
                owned.anchor()
            except FileNotFoundError:
                return
            parent = owned.parent_fd
            try:
                try:
                    fd = owned._resources.open(3, path.name, os.O_RDWR | os.O_NOFOLLOW,
                                               dir_fd=parent)
                except PermissionError:
                    fd = owned._resources.open(3, path.name, os.O_RDONLY | os.O_NOFOLLOW,
                                               dir_fd=parent)
            except FileNotFoundError:
                return
            info = os.fstat(fd)
            if not stat.S_ISREG(info.st_mode):
                raise EEDFError('record ownership is ambiguous: ' + str(path))
            try:
                record = json.loads(owned.stream.read())
            except (ValueError, UnicodeError) as error:
                raise EEDFError('record ownership is ambiguous: ' + str(path)) from error
            if (not isinstance(record, dict) or record.get('status') not in ('PASS', 'FAIL')
                    or not isinstance(record.get('run_id'), str)):
                raise EEDFError('record ownership is ambiguous: ' + str(path))
            owned.identity = (info.st_dev, info.st_ino)
            owned.run_id = record['run_id']
            owned = self._file_registry.register_record(
                str(path), owned, owned._resources, adopt=True)
        owned.revoke(payload)

    def _invalidate_acceptance(self, record_path, record):
        """Invalidate every prior authoritative record before an attempt."""
        candidates = []
        if self._acceptance_manifest_path is not None:
            candidates.append(Path(self._acceptance_manifest_path))
        if record_path is not None and Path(record_path) not in candidates:
            candidates.append(Path(record_path))
        invalid = dict(record)
        invalid.update(
            status='FAIL', comparison_passed=False,
            error='prior acceptance invalidated before qualification attempt')
        payload = self._json_bytes(invalid)
        for candidate in candidates:
            self._invalidate_record(candidate, payload)
        self._acceptance_manifest_path = None
        self._qualification_commit = None

    def _enter_failed(self, error, *, y=None, record=None,
                      record_path=None, run_id=None):
        """Make FAILED and a standing PASS mutually exclusive.

        This is the only transition into FAILED.  It invalidates the caller's
        named record, the active attempt path, and every record this provider
        has returned before attempting the diagnostic/atomic FAIL publication.
        Publication faults are retained as evidence and never undo the latch.
        """
        returned_path = self._qualification_record_path
        accepted_path = self._acceptance_manifest_path
        attempt_path = (self._qualification_attempt.get('record_path')
                        if self._qualification_attempt is not None else None)
        caller_path = None
        caller_named = isinstance(y, Mapping) and 'record_path' in y
        if caller_named:
            caller_path = y.get('record_path')

        if record is None:
            terminal_state = {}
            if caller_named:
                terminal_state['record_path'] = caller_path
            if isinstance(y, Mapping) and 'scratch_root' in y:
                terminal_state['scratch_root'] = y.get('scratch_root')
            record = {
                'status': 'FAIL',
                'run_id': run_id,
                'comparison_passed': False,
                'terminal_state': terminal_state,
            }
        else:
            record = dict(record)
        record['status'] = 'FAIL'
        record['comparison_passed'] = False
        try:
            record['error'] = str(error)
        except BaseException:
            record['error'] = type(error).__name__
        if run_id is None:
            run_id = record.get('run_id') or 'failed-{0:x}'.format(id(record))
        record['run_id'] = run_id

        raw_candidates = [
            caller_path, record_path, attempt_path, returned_path, accepted_path,
        ]
        # The first finalised failure wins. setdefault does not split the
        # decision into a Python check followed by a stale assignment after a
        # returning handler; nested failures use this same record/run identity.
        completed = False
        try:
            record = self._qualification_failure.setdefault('record', record)
            errors = self._revoke_publications()
            self._sync_qualification_responses(record)
            result = self._finalize_failure_record(
                record, raw_candidates, caller_named, errors)
            self._sync_qualification_responses(record)
            completed = True
            return result
        finally:
            if not completed:
                # A second signal can interrupt revocation or a response
                # update. Finish both idempotent operations, and the FAIL file,
                # before that interruption can return to the outer publisher.
                record = self._qualification_failure.setdefault('record', record)
                errors = self._revoke_publications()
                self._sync_qualification_responses(record)
                self._finalize_failure_record(
                    record, raw_candidates, caller_named, errors)
                self._sync_qualification_responses(record)

    def _finalize_failure_record(self, record, raw_candidates, caller_named,
                                 publication_errors):
        run_id = record['run_id']
        candidates = []
        for candidate in raw_candidates:
            if candidate is None:
                continue
            try:
                candidate = Path(candidate)
            except BaseException as path_error:
                publication_errors.append(path_error)
                continue
            if candidate not in candidates:
                candidates.append(candidate)

        try:
            payload = self._json_bytes(record)
        except BaseException as serialization_error:
            publication_errors.append(serialization_error)
            payload = self._json_bytes({
                'status': 'FAIL',
                'run_id': run_id,
                'comparison_passed': False,
                'error': 'qualification failure record was not serializable',
            })
        for candidate in candidates:
            try:
                self._invalidate_record(candidate, payload)
            except BaseException as invalidation_error:
                publication_errors.append(invalidation_error)

        self._acceptance_manifest_path = None
        self._qualification_commit = None
        self._qualification_record = record
        self._qualification_record_path = None
        self._qualification_run_id = run_id

        target = candidates[0] if candidates else None
        state = record.get('terminal_state')
        if target is None and not caller_named and isinstance(state, Mapping):
            try:
                if 'scratch_root' not in state:
                    saved_spec = json.loads(
                        (self._table.path / 'generation_spec.json').read_text())
                    state['scratch_root'] = str(saved_spec['scratch_root'])
                target = self._record_path(state, run_id)
            except BaseException as path_error:
                publication_errors.append(path_error)
        if target is not None:
            try:
                self._write_failure(record, target, run_id)
            except BaseException as publication_error:
                publication_errors.append(publication_error)
                self._qualification_record_path = None
        if publication_errors:
            record['publication_error'] = '; '.join(
                str(item) for item in publication_errors)
            return publication_errors[-1]
        return None

    @ownership_scope
    def _publish_acceptance(self, commit, record_path):
        """Publish once from a source that reentrant failure can revoke."""
        context = self._qualification_current
        if (context is None or context.commit is not commit
                or not self._publication_is_current(context)):
            return None
        owned = self._new_record_file(record_path, commit['run_id'], 'acceptance')
        context.file = owned
        owned._resources.retain()
        temporary = owned.source
        context.temporary = temporary
        try:
            if not self._publication_is_current(context):
                return None
            self._stage_record(owned, commit['manifest'])
            if not self._publication_is_current(context):
                return None
            self._validate_acceptance_manifest(
                temporary, commit['run_id'], commit, staged=True)
            if not self._publication_is_current(context):
                return None
            # After this check, reentry removes/invalidates this exact source.
            # No subsequent path recreates it or retries the stale publication.
            os.link(temporary, record_path, follow_symlinks=False)
            context.published = True
            context.generation['published'] = True
            owned.finish()
        finally:
            if not context.published:
                owned.revoke(
                    self._json_bytes({'status': 'FAIL', 'run_id': commit['run_id']}))
        return record_path

    def _validate_acceptance_manifest(self, path, run_id, commit=None, staged=False):
        """Validate files against process-local committed bytes, never manifest claims."""
        if path is None or not Path(path).is_file():
            raise EEDFError('current-run acceptance manifest is absent')
        path = Path(path)
        commit = self._qualification_commit if commit is None else commit
        if commit is None or commit.get('run_id') != run_id:
            raise EEDFError('acceptance manifest is not for the current run')
        record_path = Path(commit['record_path'])
        if not staged and record_path.resolve() != path.resolve():
            raise EEDFError('acceptance manifest was read from the wrong output location')
        artifacts = commit['artifacts']
        if tuple(artifacts) != MANDATORY_QUALIFICATION_ARTIFACTS:
            raise EEDFError('committed acceptance has an invalid artifact set')
        run_directory = Path(commit['run_directory']).resolve()
        for name in MANDATORY_QUALIFICATION_ARTIFACTS:
            candidate = (run_directory / name).resolve()
            if run_directory not in candidate.parents or not candidate.is_file():
                raise EEDFError('acceptance artifact is absent: ' + name)
            if candidate.read_bytes() != artifacts[name]:
                raise EEDFError(name + ' artifact bytes differ from committed record')
        if path.read_bytes() != commit['manifest_bytes']:
            raise EEDFError('acceptance manifest bytes differ from committed record')
        return copy.deepcopy(commit['manifest'])

    def _write_development_record(self, record, state, run_id):
        path = self._record_path(state, run_id)
        self._qualification_record = record
        self._qualification_record_path = self._write_json(path, record)
        return path

    @staticmethod
    def _tolerance_check(checks, rule, quantity, actual, reference, tolerance,
                         labels=None, relative_scale=None):
        actual = np.asarray(actual, dtype=float)
        reference = np.asarray(reference, dtype=float)
        if actual.shape != reference.shape:
            raise EEDFError('{0} shape mismatch for {1}'.format(rule, quantity))
        error = np.abs(actual - reference)
        scale = np.abs(reference) if relative_scale is None else np.asarray(
            relative_scale, dtype=float)
        allowed = (np.asarray(tolerance['atol'], dtype=float) +
                   float(tolerance['rtol']) * scale)
        allowed = np.broadcast_to(allowed, error.shape)
        passed = np.isfinite(error) & (error <= allowed)
        if error.shape:
            scores = np.full(error.shape, np.inf, dtype=float)
            finite = np.isfinite(error) & np.isfinite(allowed)
            scores[finite] = error[finite] - allowed[finite]
            flat = int(np.argmax(scores))
            index = list(np.unravel_index(flat, error.shape))
            label = labels[flat] if labels is not None and flat < len(labels) else index
            absolute_error = float(error.flat[flat])
            allowed_error = float(allowed.flat[flat])
            relative_denominator = float(np.broadcast_to(scale, error.shape).flat[flat])
        else:
            index = []
            label = labels[0] if labels else quantity
            absolute_error = float(error)
            allowed_error = float(allowed)
            relative_denominator = float(scale)
        checks.append({
            'rule': rule,
            'quantity': quantity,
            'argmax': label,
            'index': index,
            'absolute_error': absolute_error,
            'relative_error': (absolute_error / abs(relative_denominator)
                               if relative_denominator and
                               np.isfinite(absolute_error) and
                               np.isfinite(relative_denominator) else None),
            'allowed_error': allowed_error,
            'passed': bool(np.all(passed)),
        })

    @staticmethod
    def _quantity_drift(quantity, actual, reference, labels=None):
        actual = np.asarray(actual, dtype=float)
        reference = np.asarray(reference, dtype=float)
        if actual.shape != reference.shape:
            return {'quantity': quantity, 'shape_mismatch': [list(actual.shape),
                                                              list(reference.shape)]}
        error = np.abs(actual - reference)
        if error.shape:
            flat = int(np.argmax(error))
            index = list(np.unravel_index(flat, error.shape))
            argmax = labels[flat] if labels is not None and flat < len(labels) else index
            absolute = float(error.flat[flat])
            ref = float(reference.flat[flat])
        else:
            index, argmax = [], quantity
            absolute, ref = float(error), float(reference)
        return {'quantity': quantity, 'argmax': argmax, 'index': index,
                'absolute_error': absolute,
                'relative_error': (absolute / abs(ref)
                                   if ref and np.isfinite(absolute) and
                                   np.isfinite(ref) else None)}

    def _all_quantity_drifts(self, predicted, direct):
        channels = [entry.get('description', str(i))
                    for i, entry in enumerate(self._runtime_metadata['channel_map'])]
        drifts = [self._quantity_drift('EN_Td', predicted['EN_Td'], direct['EN_Td'])]
        for name in sorted(set(predicted['swarm']) | set(direct['swarm'])):
            if name in predicted['swarm'] and name in direct['swarm']:
                drifts.append(self._quantity_drift(
                    'swarm.' + name, predicted['swarm'][name], direct['swarm'][name]))
        for name in ('k_ine', 'k_sup', 'channel_power', 'attachment_energy_eV',
                     'target_fractions', 'product_fractions'):
            if name in predicted and name in direct:
                drifts.append(self._quantity_drift(name, predicted[name], direct[name], channels))
        for name in sorted(set(predicted['power_groups']) | set(direct['power_groups'])):
            if name in predicted['power_groups'] and name in direct['power_groups']:
                drifts.append(self._quantity_drift(
                    'power.' + name, predicted['power_groups'][name],
                    direct['power_groups'][name]))
        if 'f0' in predicted and 'f0' in direct:
            drifts.append(self._quantity_drift('f0', predicted['f0'], direct['f0']))
        return drifts

    def _verify_generation_spec(self, spec):
        """Verify that a saved direct-solver request owns this fingerprint."""
        manifest = self._json_value(self._runtime_metadata)
        row_inputs = manifest['row_inputs']
        if _content_hash(row_inputs) != manifest['fingerprint']:
            raise EEDFError('generation spec mismatch: manifest fingerprint')
        stored_setup = row_inputs.get('qualification_setup_sha256')
        if not stored_setup:
            raise EEDFError(
                'generation spec mismatch: artifact has no physical setup fingerprint')
        if manifest.get('qualification_setup_sha256') != stored_setup:
            raise EEDFError(
                'generation spec mismatch: physical setup fingerprint record')
        solver = manifest.get('solver')
        if not isinstance(solver, dict):
            raise EEDFError('generation spec mismatch: solver identity')

        def require_equal(label, actual, expected):
            if _canonical_json(actual) != _canonical_json(expected):
                raise EEDFError('generation spec mismatch: ' + label)

        require_equal('solver commit', spec.get('loki_commit'),
                      solver.get('commit'))
        require_equal('binary', spec.get('binary', {}).get('sha256'),
                      solver.get('binary_sha256'))
        require_equal('compiler', spec.get('compiler'), solver.get('compiler'))
        require_equal('cmake cache', spec.get('cmake_cache', {}).get('sha256'),
                      solver.get('cmake_cache_sha256'))
        spec_options = spec.get('solver_options')
        saved_options = solver.get('options')
        if (isinstance(spec_options, dict) and isinstance(saved_options, dict) and
                spec_options.get('ionizationOperatorType') !=
                saved_options.get('ionizationOperatorType')):
            raise EEDFError('generation spec mismatch: operator')
        require_equal('solver options', spec_options, saved_options)
        require_equal('working conditions', spec.get('working_conditions'),
                      solver.get('working_conditions'))
        spec_inputs = {
            name: {'sha256': item.get('sha256'), 'kind': item.get('kind')}
            for name, item in spec.get('input_files', {}).items()
        }
        require_equal('input SHAs', spec_inputs, solver.get('input_files'))

        cross_sections = {name: item['sha256']
                          for name, item in spec_inputs.items()
                          if item['kind'] == 'cross_section'}
        auxiliary = {name: item['sha256']
                     for name, item in spec_inputs.items()
                     if item['kind'] != 'cross_section'}
        for label, key, value in (
                ('solver commit', 'loki_commit', spec.get('loki_commit')),
                ('binary', 'binary_sha256',
                 spec.get('binary', {}).get('sha256')),
                ('solver options', 'solver_options', spec_options),
                ('working conditions', 'working_conditions',
                 spec.get('working_conditions')),
                ('qualification quantity rules', 'qualification_quantity_rules',
                 spec.get('qualification_quantity_rules')),
                ('input SHAs', 'cross_sections', cross_sections),
                ('input SHAs', 'auxiliary_inputs', auxiliary)):
            if key not in row_inputs:
                raise EEDFError(
                    'generation spec mismatch: fingerprint has no ' + key)
            require_equal(label, value, row_inputs[key])
        return stored_setup

    @staticmethod
    def _declared_values(declarations, label):
        """Parse setup declarations numerically while retaining their names."""
        if not isinstance(declarations, list):
            raise EEDFError(label + ' declarations are invalid')
        values = {}
        pattern = re.compile(
            r'^\s*(.+?)\s*=\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*$')
        for declaration in declarations:
            match = pattern.fullmatch(str(declaration))
            if match is None:
                raise EEDFError(label + ' declaration is invalid')
            name, value = match.groups()
            value = float(value)
            if name in values or not np.isfinite(value):
                raise EEDFError(label + ' declaration is invalid')
            values[name] = value
        return values

    @staticmethod
    def _terminal_values(state, key, expected):
        values = state.get(key)
        if not isinstance(values, Mapping) or set(values) != set(expected):
            raise EEDFError('terminal ' + key.replace('_', ' ') + ' are incomplete')
        result = {str(name): float(value) for name, value in values.items()}
        if not all(np.isfinite(value) for value in result.values()):
            raise EEDFError('terminal ' + key.replace('_', ' ') + ' are incomplete')
        return result

    @staticmethod
    def _validate_qualification_row(row, label):
        """Reject an incomplete, non-finite, or non-converged solver row."""
        required = {
            'u', 'EN_Td', 'composition', 'swarm', 'gas_fractions',
            'state_populations', 'k_ine', 'k_sup', 'channel_power',
            'power_groups', 'attachment_energy_eV', 'target_fractions',
            'product_fractions', 'rate_floors', 'below_floor', 'f0',
            'energy_eV', 'energy_edges_eV', 'converged', 'iteration_count',
        }
        if not isinstance(row, Mapping):
            raise EEDFError('invalid ' + label + ' row: not a mapping')
        missing = sorted(required - set(row))
        if missing:
            raise EEDFError(
                'invalid ' + label + ' row: missing ' + ', '.join(missing))
        if row['converged'] is not True:
            raise EEDFError('invalid ' + label + ' row: converged')
        mappings = ('composition', 'swarm', 'gas_fractions',
                    'state_populations', 'power_groups')
        for name in mappings:
            values = row[name]
            if not isinstance(values, Mapping) or not values:
                if name == 'composition' and isinstance(values, Mapping):
                    continue
                raise EEDFError('invalid ' + label + ' row: ' + name)
            try:
                finite = all(np.all(np.isfinite(np.asarray(value, dtype=float)))
                             for value in values.values())
            except (TypeError, ValueError):
                finite = False
            if not finite:
                raise EEDFError('invalid ' + label + ' row: ' + name)
        numeric = required - set(mappings) - {'converged'}
        for name in numeric:
            try:
                finite = np.all(np.isfinite(np.asarray(row[name], dtype=float)))
            except (TypeError, ValueError):
                finite = False
            if not finite:
                raise EEDFError('invalid ' + label + ' row: ' + name)

    def _check_all_recorded_quantities(self, checks, predicted, direct, tolerances):
        """Apply only a-posteriori A budgets to every supplementary row value."""
        channels = [entry.get('description', str(i))
                    for i, entry in enumerate(self._runtime_metadata['channel_map'])]
        rules = self._runtime_metadata['row_inputs'].get(
            'qualification_quantity_rules')
        if not isinstance(rules, Mapping):
            raise EEDFError('artifact has no qualification quantity rules')

        def check(quantity, actual, reference, labels=None):
            rule = rules.get(quantity)
            allowed = qualification_rule_for_quantity(quantity)
            if rule != allowed:
                if rule is not None and allowed is not None:
                    raise EEDFError(
                        'invalid qualification rule for ' + quantity +
                        ': expected ' + allowed + ', got ' + str(rule))
                raise EEDFError(
                    'artifact has no qualification rule for ' + quantity)
            self._tolerance_check(checks, rule, quantity, actual, reference,
                                  tolerances[rule], labels)

        check('EN_Td', predicted['EN_Td'], direct['EN_Td'])
        predicted_names = set(predicted['swarm'])
        direct_names = set(direct['swarm'])
        if predicted_names != direct_names:
            missing = sorted(predicted_names ^ direct_names)
            raise EEDFError(
                'A2 quantity set mismatch for swarm: ' + ', '.join(missing))
        for name in sorted(predicted_names - {
                'mean_energy_eV', 'mobility_N', 'diffusion_N'}):
            check('swarm.' + name, predicted['swarm'][name],
                  direct['swarm'][name])
        predicted_names = set(predicted['power_groups'])
        direct_names = set(direct['power_groups'])
        if predicted_names != direct_names:
            missing = sorted(predicted_names ^ direct_names)
            raise EEDFError(
                'A4 quantity set mismatch for power_groups: ' + ', '.join(missing))
        for name in sorted(predicted_names):
            check('power.' + name, predicted['power_groups'][name],
                  direct['power_groups'][name])
        if 'channel_power' not in predicted or 'channel_power' not in direct:
            raise EEDFError('A5 quantity set mismatch: channel_power')
        predicted_field = abs(float(predicted['power_groups']['field']))
        direct_field = abs(float(direct['power_groups']['field']))
        if not predicted_field or not direct_field:
            raise EEDFError('A5 cannot scale zero field power')
        check('channel_power_fraction',
              np.asarray(predicted['channel_power']) / predicted_field,
              np.asarray(direct['channel_power']) / direct_field, channels)
        for name in ('attachment_energy_eV', 'target_fractions',
                     'product_fractions'):
            if name not in predicted or name not in direct:
                raise EEDFError('qualification quantity set mismatch: ' + name)
            check(name, predicted[name], direct[name], channels)
        if any(name not in predicted or name not in direct
               for name in ('f0', 'energy_eV', 'energy_edges_eV')):
            raise EEDFError('A1 quantity set mismatch: f0')
        if (not np.array_equal(predicted['energy_eV'], direct['energy_eV']) or
                not np.array_equal(predicted['energy_edges_eV'],
                                   direct['energy_edges_eV'])):
            raise EEDFError('A5 f0 energy grid mismatch')
        weights = (np.sqrt(np.asarray(direct['energy_eV'], dtype=float)) *
                   np.diff(np.asarray(direct['energy_edges_eV'], dtype=float)))
        distance = np.sum(
            np.abs(np.asarray(predicted['f0'], dtype=float) -
                   np.asarray(direct['f0'], dtype=float)) * weights)
        check('f0.weighted_L1', distance, 0.)

    def _compare_qualification_rows(self, predicted, direct, state, tolerances):
        required = tuple('A' + str(i) for i in range(1, 7))
        missing = [name for name in required if name not in tolerances]
        if missing:
            raise EEDFError('artifact has no runtime tolerance ' + ', '.join(missing))
        checks = []
        self._check_all_recorded_quantities(
            checks, predicted, direct, tolerances)
        channels = [entry.get('description', str(i))
                    for i, entry in enumerate(self._runtime_metadata['channel_map'])]
        self._tolerance_check(checks, 'A1', 'mean_energy_eV',
                              predicted['swarm']['mean_energy_eV'],
                              direct['swarm']['mean_energy_eV'], tolerances['A1'])
        swarm_names = ('mobility_N', 'diffusion_N')
        for name in swarm_names:
            if name not in predicted['swarm'] or name not in direct['swarm']:
                raise EEDFError('A2 quantity set mismatch: ' + name)
            self._tolerance_check(checks, 'A2', 'swarm.' + name,
                                  predicted['swarm'][name], direct['swarm'][name],
                                  tolerances['A2'])
        for name in ('k_ine', 'k_sup'):
            identities = []
            for index, channel in enumerate(self._runtime_metadata['channel_map']):
                match = re.fullmatch(
                    r'e \+ (.+?) (?:<->|->) e \+ (.+), [A-Za-z]+',
                    channel['description'])
                if match:
                    source, product = match.groups()
                else:
                    source = product = channel.get(
                        'flux_group', channel['description'])
                if name == 'k_sup':
                    source, product = product, source
                identities.append((source, product))
            floors = self._runtime_metadata['floors']

            def materiality(row):
                fractions = np.asarray(
                    row['target_fractions'] if name == 'k_ine'
                    else row['product_fractions'], dtype=float)
                rates = np.asarray(row[name], dtype=float)
                fluxes = np.abs(fractions * rates)
                losses = {}
                productions = {}
                for index, (source, product) in enumerate(identities):
                    losses[source] = losses.get(source, 0.) + fluxes[index]
                    productions[product] = (
                        productions.get(product, 0.) + fluxes[index])
                rate_floors = np.asarray(predicted['rate_floors'], dtype=float)
                row_floors = np.asarray(row['rate_floors'], dtype=float)
                if (rate_floors.shape != row_floors.shape or
                        not np.array_equal(rate_floors, row_floors)):
                    raise EEDFError('A3 physical rate floors differ from artifact')
                result = []
                for index, (source, product) in enumerate(identities):
                    totals = [value for value in
                              (losses[source], productions[product]) if value > 0.]
                    scale = min(totals, default=0.)
                    flux_floor = (
                        float(floors['absolute_flux_fraction']) * scale /
                        abs(fractions[index]) if fractions[index] else 0.)
                    load_bearing = (
                        scale > 0. and
                        fluxes[index] >=
                        float(floors['rate_flux_fraction']) * scale and
                        abs(rates[index]) >= rate_floors[index])
                    result.append((load_bearing, rate_floors[index], flux_floor))
                return result

            predicted_materiality = materiality(predicted)
            direct_materiality = materiality(direct)
            for index, channel in enumerate(self._runtime_metadata['channel_map']):
                classification = channel.get('classification')
                predicted_load = predicted_materiality[index][0]
                direct_load = direct_materiality[index][0]
                if classification != 'A' and not (
                        classification == 'B' and
                        (predicted_load or direct_load)):
                    continue
                tolerance = dict(tolerances['A3'])
                tolerance['atol'] = max(
                    float(tolerance['atol']),
                    float(predicted_materiality[index][1]),
                    float(direct_materiality[index][1]),
                    float(predicted_materiality[index][2]),
                    float(direct_materiality[index][2]))
                if not predicted_load and not direct_load:
                    tolerance['rtol'] = 0.
                self._tolerance_check(
                    checks, 'A3', name + ':' + str(index),
                    predicted[name][index], direct[name][index], tolerance,
                    [channels[index]])
        loss_names = ('elastic_gain', 'elastic_loss', 'car_gain', 'car_loss',
                      'excitation_loss', 'excitation_gain', 'vibrational_loss',
                      'vibrational_gain', 'rotational_loss', 'rotational_gain',
                      'ionization', 'attachment')
        self._tolerance_check(
            checks, 'A4', 'power.total_loss',
            -sum(float(predicted['power_groups'].get(name, 0.)) for name in loss_names),
            -sum(float(direct['power_groups'].get(name, 0.)) for name in loss_names),
            tolerances['A4'])
        field = abs(float(direct['power_groups']['field']))
        a5 = tolerances['A5']
        significant = [i for i, value in enumerate(direct['channel_power'])
                       if abs(value) >= float(a5.get('share_min', 0.)) * field]
        for index in significant:
            predicted_fraction = (float(predicted['channel_power'][index]) /
                                  abs(float(predicted['power_groups']['field'])))
            direct_fraction = float(direct['channel_power'][index]) / field
            self._tolerance_check(checks, 'A5', 'channel_power_fraction:' + str(index),
                                  predicted_fraction, direct_fraction, a5,
                                  [channels[index]])
        budget = state.get('energy_budget', {})
        if budget:
            predicted_field = float(predicted['power_groups']['field'])
            if predicted_field == 0.:
                raise EEDFError('A6 cannot scale a zero interpolated field power')
            scale = float(budget['joule_power']) / predicted_field
            direct_joule = scale * float(direct['power_groups']['field'])
            engine_target = (float(budget['P_abs']) -
                             float(budget.get('Q_wall_electron', 0.)) -
                             float(budget.get('Q_wall_ion', 0.)) -
                             float(budget.get('Q_flow', 0.)) +
                             scale * float(direct['power_groups'].get('growth', 0.)))
            self._tolerance_check(checks, 'A6', 'joule_power',
                                  direct_joule, engine_target,
                                  tolerances['A6'],
                                  relative_scale=abs(float(budget['P_abs'])))
        else:
            raise EEDFError('terminal qualification state has no energy_budget')
        return checks

    @_qualification_locked
    def _run_qualification(self, y, t, runner, development):
        from rmgpy.exceptions import PlasmaStateError
        from rmgpy.tools.eedf.integrity import _resolved_properties
        from rmgpy.tools.eedf.loki import LoKIDriver, _enrich_row

        start = None
        run_id = None
        record_path = (None if development else
                       self._acceptance_manifest_path)
        state = {}
        record = {
            'status': DEVELOPMENT_DIAGNOSTIC_STATUS if development else 'FAIL',
            'run_id': None,
            'comparison_passed': False,
            'terminal_state': state,
        }
        if not development:
            self._qualification_attempt = {
                'record': record,
                'record_path': record_path,
            }
        try:
            if not development:
                # This is deliberately the first fallible action in a terminal
                # attempt.  A prior PASS must stop being authoritative before
                # clocks, identifiers, state callbacks, or paths are touched.
                self._invalidate_acceptance(None, record)
                self._qualification_phase = None
                self._qualification_record = record
                self._qualification_record_path = None
            self._acceptance_manifest_path = None
            self._qualification_commit = None
            start = time.monotonic()
            run_id = uuid.uuid4().hex
            record['run_id'] = run_id
            self._qualification_run_id = run_id
            record.update(
                artifact_sha256=self._runtime_metadata['artifact_sha256'],
                generation_spec=str(
                    self._table.path / 'generation_spec.json'),
                branch=self._branch,
            )
            if isinstance(y, dict) and 'record_path' in y:
                state['record_path'] = y['record_path']
                if not development:
                    self._capture_qualification_state(state)
            state['t'] = float(t)
            state = self._qualification_state(y, t)
            record['terminal_state'] = state
            if not development:
                record_path = self._qualification_attempt['record_path']
            if state.get('record_path') is None and 'record_path' in state:
                raise EEDFError(
                    'terminal qualification state has invalid record_path')
            if 'record_path' in state:
                record_path = Path(state['record_path'])
            predicted = self.row(state['u'], state['composition'])._as_dict()
            self._validate_qualification_row(predicted, 'interpolated')
            record['interpolated_row'] = predicted
            spec_path = self._table.path / 'generation_spec.json'
            spec = json.loads(spec_path.read_text())
            stored_table_identity = self._verify_generation_spec(spec)
            if 'scratch_root' in state:
                spec['scratch_root'] = str(state['scratch_root'])
            if 'scratch_root' not in spec:
                raise EEDFError('qualification execution has no scratch_root')
            state['scratch_root'] = str(spec['scratch_root'])
            record_path = self._record_path(state, run_id)
            gas_properties, state_properties, saved_gas, saved_populations = (
                _resolved_properties(spec, state['composition']))
            terminal_gas = self._terminal_values(
                state, 'gas_fractions', saved_gas)
            terminal_populations = self._terminal_values(
                state, 'state_populations', saved_populations)
            gas_properties['fraction'] = [
                '{0} = {1:.17g}'.format(name, terminal_gas[name])
                for name in sorted(terminal_gas)]
            state_properties['population'] = [
                '{0} = {1:.17g}'.format(name, terminal_populations[name])
                for name in sorted(terminal_populations)]
            if 'statisticalWeight' in state_properties:
                saved_weights = self._declared_values(
                    state_properties['statisticalWeight'],
                    'state statistical weights')
                terminal_weights = self._terminal_values(
                    state, 'state_statistical_weights', saved_weights)
                if any(terminal_weights[name] != saved_weights[name]
                       for name in saved_weights):
                    raise EEDFError(
                        'generation spec mismatch: physical setup fingerprint '
                        '(terminal statistical weights)')
            spec['gas_properties'] = gas_properties
            spec['state_properties'] = state_properties
            terminal_fields = [float(np.exp(state['u']))]
            job_name = 'qualification-' + run_id
            if runner is None:
                binary = spec.get('binary', {}).get('path')
                if not binary or not Path(binary).is_file() or not os.access(binary, os.X_OK):
                    raise EEDFError('missing or unexecutable LoKI-B binary')
                driver = LoKIDriver(spec)
            else:
                driver = object.__new__(LoKIDriver)
                driver.spec = spec
                driver.setups = {}
            bundle = driver.prepare(
                state['composition'], terminal_fields, job_name,
                terminal_state=state, verify_dependencies=runner is None)
            if bundle.table_identity != stored_table_identity:
                raise EEDFError(
                    'generation spec mismatch: physical setup fingerprint')
            if runner is None:
                raw = driver.run_prepared(bundle)[0]
                channel_map = self._json_value(self._runtime_metadata['channel_map'])
                direct = _enrich_row(raw, channel_map, bundle.spec)
            else:
                direct = runner(bundle)[0]
                work = Path(bundle.spec['scratch_root']) / bundle.job_name
                work.mkdir(parents=True, exist_ok=False)
                (work / 'setup.in').write_bytes(bundle.setup_bytes)
            record['solver_row'] = direct
            record['table_identity'] = bundle.table_identity
            record['execution_identity'] = bundle.execution_identity
            terminal_identity = _content_hash(self._json_value(state))
            record['terminal_state_identity'] = terminal_identity
            if direct.get('execution_identity') != bundle.execution_identity:
                raise EEDFError(
                    'direct solver output execution identity differs from prepared bundle')
            self._validate_qualification_row(direct, 'direct')
            tolerances = self._runtime_metadata['tolerances']
            if development and not all('A' + str(i) in tolerances for i in range(1, 7)):
                # Historical development artifacts predate A1--A6. Reuse their
                # explicit held-out budgets for measurement only; this route can
                # never confer qualification.
                tolerances = dict(tolerances)
                for target, source in (('A1', 'H1'), ('A2', 'H2'), ('A3', 'H3'),
                                       ('A4', 'H4'), ('A5', 'H4')):
                    if target not in tolerances and source in tolerances:
                        tolerances[target] = dict(tolerances[source])
                if 'A5' in tolerances:
                    tolerances['A5'].setdefault(
                        'share_min', self._runtime_metadata.get('floors', {}).get(
                            'relative_power_share', 0.))
                if 'A6' not in tolerances and 'energy_budget' in state:
                    budget = state['energy_budget']
                    tolerances['A6'] = {'rtol': float(budget['A6b-runtime']['tolerance']
                                                    if 'A6b-runtime' in budget else
                                                    budget['A6b_tolerance']), 'atol': 0.}
                record['tolerance_source'] = 'development artifact held-out budgets'
            checks = self._compare_qualification_rows(predicted, direct, state, tolerances)
            record.update(checks=checks,
                          quantity_drifts=self._all_quantity_drifts(predicted, direct),
                          drift_maxima={rule: max(
                              (item for item in checks if item['rule'] == rule),
                              key=lambda item: (item['absolute_error']
                                                if np.isfinite(item['absolute_error'])
                                                else float('inf')),
                              default=None)
                              for rule in sorted({item['rule'] for item in checks})})
            budget = state.get('energy_budget', {})
            if 'A6b-runtime' in budget:
                record['engine_A6b'] = self._json_value(budget['A6b-runtime'])
            elif 'A6b_relative' in budget:
                record['engine_A6b'] = {
                    'relative_error': float(budget['A6b_relative']),
                    'tolerance': float(budget['A6b_tolerance']),
                    'passed': bool(float(budget['A6b_relative']) <=
                                   float(budget['A6b_tolerance'])),
                }
            record['comparison_passed'] = all(item['passed'] for item in checks)
            record['wall_time_s'] = time.monotonic() - start
            if development:
                self.annotate_development_diagnostic(record)
                self._write_development_record(record, state, run_id)
                return record
            failed = next((item for item in checks if not item['passed']), None)
            if failed is not None:
                raise PlasmaStateError('EEDF qualification failed: ' + failed['rule'])
            work = Path(bundle.spec['scratch_root']) / bundle.job_name
            record['status'] = 'CANDIDATE'
            artifact_bytes = {
                'setup.in': bundle.setup_bytes,
                'interpolated-row.json': self._json_bytes(predicted),
                'solver-row.json': self._json_bytes(direct),
            }
            self._write_json(work / 'interpolated-row.json', predicted)
            self._write_json(work / 'solver-row.json', direct)
            self._write_json(work / 'qualification-evidence.json', record)
            artifact_bytes['qualification-evidence.json'] = self._json_bytes(record)
            artifacts = {
                name: hashlib.sha256(artifact_bytes[name]).hexdigest()
                for name in MANDATORY_QUALIFICATION_ARTIFACTS
            }
            manifest = dict(record)
            manifest.update(
                status='PASS', artifacts=artifacts,
                run_directory=str(work.resolve()),
                acceptance_manifest=str(record_path.resolve()))
            committed_manifest = self._json_value(manifest)
            commit = MappingProxyType({
                'run_id': run_id,
                'record_path': str(record_path.resolve()),
                'run_directory': str(work.resolve()),
                'artifacts': MappingProxyType({
                    name: artifact_bytes[name]
                    for name in MANDATORY_QUALIFICATION_ARTIFACTS
                }),
                'manifest': committed_manifest,
                'manifest_bytes': self._json_bytes(committed_manifest),
            })
            context = self._qualification_current
            context.commit = commit
            context.record_path = record_path
            self._sync_qualification_responses(manifest)
            self._qualification_record_path = record_path
            self._acceptance_manifest_path = record_path
            self._qualification_commit = commit
            # Keep the candidate in memory through signal release. The outer
            # attempt publishes only after every owned cancellation has been
            # accounted for; no PASS file can survive death during restoration.
            self._qualification_phase = 'CANDIDATE'
            return context.record
        except BaseException as error:
            if start is not None:
                try:
                    record['wall_time_s'] = time.monotonic() - start
                except BaseException:
                    pass
            record['error'] = str(error)
            record['comparison_passed'] = False
            if development:
                self.annotate_development_diagnostic(record)
                if record_path is None:
                    record_path = self._record_path(state, run_id)
                self._write_development_record(record, state, run_id)
                return record
            if self._qualification_attempt is not None:
                record_path = self._qualification_attempt['record_path']
            publication_error = self._enter_failed(
                error, y=y, record=record, record_path=record_path,
                run_id=run_id)
            if isinstance(error, PlasmaStateError):
                raise
            if not isinstance(error, Exception):
                raise
            if (publication_error is not None
                    and not isinstance(publication_error, Exception)):
                raise publication_error
            if publication_error is not None:
                error = publication_error
            raise PlasmaStateError('EEDF qualification failed: ' + str(error)) from error
        finally:
            if not development:
                self._qualification_attempt = None

    @_qualification_attempt_locked
    def record_terminal_failure(self, y, t, error):
        """Replace any prior PASS with one recorded, terminal FAILED latch."""
        run_id = None
        record_path = self._acceptance_manifest_path
        state = {}
        record = {
            'status': 'FAIL',
            'run_id': None,
            'comparison_passed': False,
            'terminal_state': state,
        }
        self._qualification_attempt = {
            'record': record,
            'record_path': record_path,
        }
        pending = None
        try:
            # Keep invalidation, state resolution, FAILED publication, and the
            # latch under this one lock.  Interruptions are re-raised only once
            # the attempt is durably non-PASS and the in-memory latch is final.
            self._invalidate_acceptance(None, record)
            self._qualification_record = record
            self._qualification_record_path = None
            run_id = uuid.uuid4().hex
            record['run_id'] = run_id
            self._qualification_run_id = run_id
            record.update(
                artifact_sha256=self._runtime_metadata['artifact_sha256'],
                generation_spec=str(
                    self._table.path / 'generation_spec.json'),
                branch=self._branch,
                error=str(error),
            )
            if isinstance(y, dict) and 'record_path' in y:
                state['record_path'] = y['record_path']
                self._capture_qualification_state(state)
            state = self._qualification_state(y, t)
            record['terminal_state'] = state
        except BaseException as attempt_error:
            pending = attempt_error
            if isinstance(attempt_error, Exception):
                try:
                    state['t'] = float(t)
                except BaseException as state_error:
                    if not isinstance(state_error, Exception):
                        pending = state_error
            record['error'] = str(error)
        finally:
            if self._qualification_attempt is not None:
                record_path = self._qualification_attempt['record_path']
            publication_error = self._enter_failed(
                error, y=y, record=record, record_path=record_path,
                run_id=run_id)
            self._qualification_attempt = None
        if pending is not None and not isinstance(pending, Exception):
            raise pending
        if (publication_error is not None
                and not isinstance(publication_error, Exception)):
            raise publication_error
        return self._qualification_record

    @_qualification_attempt_locked
    def qualify(self, y, t):
        """Qualify one accepted terminal state with an independent direct solve."""
        from rmgpy.exceptions import PlasmaStateError

        if self._development_blockers:
            raise EEDFError(DEVELOPMENT_UNQUALIFIED_STATUS + ' refuses qualification')
        if self._qualification_status == 'FAILED':
            raise PlasmaStateError(
                'EEDF qualification failed: provider is latched FAILED')
        return self._run_qualification(y, t, None, development=False)

    def development_diagnostic(self, y, t, *, runner=None):
        """Measure direct/interpolated drift without changing development status."""
        if not self._development_blockers:
            raise EEDFError('development diagnostic requires a development provider')
        return self._run_qualification(y, t, runner, development=True)

    def annotate_development_diagnostic(self, diagnostic):
        """Add the mandatory development warning to one diagnostic mapping."""
        if not isinstance(diagnostic, dict):
            raise TypeError('diagnostic must be a mapping')
        notice = self.development_notice
        if notice is not None:
            diagnostic.update(notice)
        return _detached_copy(diagnostic, 'diagnostic')

    @property
    def axis_names(self):
        return tuple(self._table.axis_names)

    def axis(self, name):
        if not isinstance(name, str):
            raise TypeError('axis name must be a string')
        name = str(name)
        table = self._table
        axis_names = table.axis_names
        try:
            index = axis_names.index(name)
        except ValueError as error:
            raise KeyError(name) from error
        result = table.axes[index].copy()
        result.setflags(write=False)
        return result

    @property
    def _cache_identity_snapshot(self):
        with self._cache_lock:
            return self._cache_identity

    @property
    def _cached_row_snapshot(self):
        with self._cache_lock:
            return self._cached_row

    def _identity(self, u, composition):
        if not isinstance(composition, dict):
            raise OutOfDomain('composition must be a named coordinate mapping')
        return (float(u), tuple(sorted((name, float(value))
                                       for name, value in composition.items())))

    def _query_snapshot(self, u, composition, table=None):
        """Own and normalize one query for checking, identity and evaluation."""
        table = self._table if table is None else table
        u, composition = _detached_copy(
            (u, composition), 'EEDF query', coordinates=True)
        if not isinstance(composition, dict):
            raise OutOfDomain('composition must be a named coordinate mapping')
        supplied = tuple(composition)
        point = table._domain_check(u, composition)
        identity = {name: composition[name] for name in supplied}
        return point[0], composition, identity

    def domain_check(self, u, composition, y=None, context=None):
        """Validate an accepted state; y/context are reactor-hook context."""
        u, composition = _detached_copy(
            (u, composition), 'EEDF domain query', coordinates=True)
        table = self._table
        axis_names = tuple(table.axis_names)
        axes = tuple(table.axes)
        try:
            normalized_u, normalized_composition, _ = self._query_snapshot(
                u, composition, table)
            return table._domain_check(
                normalized_u, normalized_composition)
        except OutOfDomain:
            values = [u]
            if isinstance(composition, dict):
                values.extend(composition.get(name) for name in axis_names[1:])
            if len(values) == len(axis_names):
                for name, axis, value in zip(axis_names, axes, values):
                    if (value is not None and
                            (not np.isfinite(value) or not axis[0] <= value <= axis[-1])):
                        raise OutOfDomain(
                            '{}={!r} outside [{!r}, {!r}]'.format(
                                name, value, axis[0], axis[-1]))
            raise

    def row(self, u, composition, y=None, context=None):
        """Return this provider's immutable row handle for one full identity."""
        table = self._table
        branch = self._branch
        u, composition, identity_composition = self._query_snapshot(
            u, composition, table)
        with self._cache_lock:
            identity = self._identity(u, identity_composition)
            cached = self._cached_row
            if identity == self._cache_identity and cached is not None:
                return cached
            local_row = _freeze_row(table.row(u, composition, branch))
            if local_row.branch_id != branch:
                raise AmbiguousBranch('row branch differs from declared branch')
            self._issued_rows[id(local_row)] = local_row
            # Publish the row and its full identity as one critical section.
            self._cached_row = local_row
            self._cache_identity = identity
            return local_row

    def _require_row(self, row):
        with self._cache_lock:
            if row is None:
                row = self._cached_row
            if row is None:
                raise RuntimeError('EEDF row has not been evaluated')
            if self._issued_rows.get(id(row)) is not row:
                raise FingerprintMismatch('row was not issued by this EEDF provider')
            if row.branch_id != self._branch:
                raise AmbiguousBranch('row branch differs from declared branch')
            if row.fingerprint != self._runtime_metadata['fingerprint']:
                raise FingerprintMismatch('row fingerprint')
            return row

    def reaction_rate(self, reaction_index, row):
        """Return a mapped SI rate from an explicit provider-issued row."""
        try:
            value = np.asarray(reaction_index)
            if (value.ndim != 0 or value.dtype.kind not in 'iu' or
                    value.dtype.kind == 'b'):
                raise ValueError
            reaction_index = int(value)
        except (TypeError, ValueError, OverflowError) as error:
            raise KeyError('reaction index must be an integer scalar') from error
        if reaction_index not in self._reaction_map:
            raise KeyError('unmapped reaction index ' + repr(reaction_index))
        row = self._require_row(row)
        column, side = self._reaction_map[reaction_index]
        return float(getattr(row, 'k_' + side)[column])

    def transport(self, row):
        """Return immutable reduced transport from an explicit row handle."""
        return self._require_row(row).swarm

    def mean_energy_eV(self, row):
        """Return mean electron energy from an explicit row handle."""
        row = self._require_row(row)
        try:
            return float(row.swarm['mean_energy_eV'])
        except KeyError as error:
            raise EEDFError('row has no mean_energy_eV') from error
