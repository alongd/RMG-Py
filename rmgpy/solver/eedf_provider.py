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
from contextlib import contextmanager
from contextvars import ContextVar
import threading
from types import MappingProxyType
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
)


DEVELOPMENT_UNQUALIFIED_STATUS = 'DEVELOPMENT \u2014 UNQUALIFIED TABLE'
_development_unqualified_route = ContextVar(
    'eedf_development_unqualified_route', default=False)


class DevelopmentWallBudgetExceeded(RuntimeError):
    """The private development runner exhausted its configured wall budget."""


@contextmanager
def development_unqualified_table_route():
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


def unqualified_artifact_blockers(manifest):
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


class EEDFProvider:
    """Immutable reactor-owned view of one accepted table branch.

    Rows returned by this provider are immutable handles issued by this exact
    provider. Consumers cannot pass a row from another table or provider into
    rate, transport, or energy accessors.
    """

    def __init__(self, path, current_model, *, artifact_sha256, reaction_map,
                 branch, empirical_laws=None):
        development = bool(_development_unqualified_route.get())
        table = EEDFTable.load(path, current_model,
                               artifact_sha256=artifact_sha256,
                               require_accepted=not development)
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
        self._development_blockers = (unqualified_artifact_blockers(table.manifest)
                                      if development else ())
        self._cache_identity = None
        self._cached_row = None
        self._issued_rows = weakref.WeakValueDictionary()
        self._cache_lock = threading.RLock()

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
    def branch(self):
        return self._branch

    @property
    def reaction_map(self):
        return self._reaction_map

    @property
    def empirical_laws(self):
        return self._empirical_laws

    @property
    def manifest(self):
        """A defensive copy of the checked artifact manifest."""
        manifest = copy.deepcopy(self._table.manifest)
        notice = self.development_notice
        if notice is not None:
            manifest.update(notice)
        return manifest

    @property
    def runtime_metadata(self):
        """Immutable metadata for reactor-internal callback use."""
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

    def require_qualified(self, operation):
        """Refuse export/qualification consumers for a development provider."""
        if self._development_blockers:
            raise EEDFError(
                '{0} refuses {1}: {2}'.format(
                    DEVELOPMENT_UNQUALIFIED_STATUS, operation,
                    '; '.join(self._development_blockers)))

    def annotate_development_diagnostic(self, diagnostic):
        """Add the mandatory development warning to one diagnostic mapping."""
        notice = self.development_notice
        if notice is not None:
            diagnostic.update(notice)
        return diagnostic

    @property
    def axis_names(self):
        return tuple(self._table.axis_names)

    def axis(self, name):
        try:
            index = self._table.axis_names.index(name)
        except ValueError as error:
            raise KeyError(name) from error
        result = self._table.axes[index].copy()
        result.setflags(write=False)
        return result

    @property
    def cache_identity(self):
        with self._cache_lock:
            return self._cache_identity

    @property
    def cached_row(self):
        with self._cache_lock:
            return self._cached_row

    def _identity(self, u, composition):
        if not isinstance(composition, dict):
            raise OutOfDomain('composition must be a named coordinate mapping')
        return (float(u), tuple(sorted((name, float(value))
                                       for name, value in composition.items())))

    def domain_check(self, u, composition, y=None, context=None):
        """Validate an accepted state; y/context are reactor-hook context."""
        return self._table.domain_check(u, composition)

    def row(self, u, composition, y=None, context=None):
        """Return this provider's immutable row handle for one full identity."""
        with self._cache_lock:
            self.domain_check(u, composition, y=y, context=context)
            identity = self._identity(u, composition)
            cached = self._cached_row
            if identity == self._cache_identity and cached is not None:
                return cached
            local_row = _freeze_row(self._table.row(u, composition, self._branch))
            if local_row.branch_id != self._branch:
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
