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

"""
Quarantine of kinetics database data that must not enter a quantitative mechanism.

Some data in the database is real, provenance-bearing, and worth keeping, yet its
*evaluation* has lost the semantics of its source model -- the number the rate law
returns is not the number the data means. Deleting it destroys evidence; leaving it
in circulation puts a meaningless rate into a mechanism that then reports success.
Quarantine is the third option: keep the data, and refuse loudly at the boundary
where it would become a quantitative claim.

**The marker lives in the database, next to the data it describes.** A kinetics
family directory may carry a ``quarantine.py`` sidecar alongside ``groups.py`` and
``rules.py``::

    name = "Some_Family/quarantine"
    state = "QUARANTINED FOR QUANTITATIVE PLASMA USE"
    appliesToKineticsClass = "Marcus"
    reason = "electrochemical reference/domain unavailable"

Nothing in this module names a family, an entry, or a rate. Which families are
quarantined and which of their entries are affected are both answers that live in
the database and are computed from it.

**The manifest declares a criterion, never a list of entries.** The affected set is
computed by evaluating ``appliesToKineticsClass`` against whatever the database
currently holds, so it cannot fall out of step with the data: adding an entry
quarantines it automatically, refitting one to a different kinetics class releases
it automatically, and renumbering or relabelling entries changes nothing. A
hard-coded list of affected entries -- in the code or in the manifest -- is exactly
the failure mode this avoids.

**Scope.** The quarantine is a property of a *family*, so it does not reach ordinary
chemistry: a database with no ``quarantine.py`` anywhere behaves precisely as it did
before, and the gate costs one ``None`` check per reaction.

**What is deliberately *not* gated.** The database still loads, the family stays
registered, ``KineticsFamily.generate_reactions`` still returns reactions, and the
rate law itself -- ``Marcus.get_rate_coefficient``, ``Reaction.get_rate_coefficient``
-- is untouched. Only admission to an RMG reaction model refuses. A consumer that
supplies the missing reference (an electrode/electrolyte model, a refit, an Arkane
job, a provenance audit) keeps working on the same unmodified data.
"""

import logging
import os.path

import rmgpy.kinetics
from rmgpy.exceptions import DatabaseError, QuarantinedKineticsError

#: Filename of the sidecar manifest inside a kinetics family directory.
QUARANTINE_FILENAME = 'quarantine.py'

#: Manifest fields without which a quarantine cannot say what it covers or why.
#: ``name``, ``shortDesc`` and ``longDesc`` are optional and free text; they follow the
#: rest of the database in being camelCase in the data file and snake_case on the object.
_REQUIRED_FIELDS = ('state', 'appliesToKineticsClass', 'reason')


class KineticsQuarantine(object):
    """
    A machine-readable record that a kinetics family's data must not enter a
    quantitative mechanism.

    ============================ ===========================================================
    Attribute                    Description
    ============================ ===========================================================
    `family_label`               Label of the family the manifest was found in
    `state`                      The campaign state string, e.g. ``"QUARANTINED FOR ..."``
    `kinetics_class`             The kinetics *class* the quarantine applies to
    `kinetics_class_name`        The name as written in the manifest
    `reason`                     Why the data cannot be evaluated meaningfully here
    `path`                       Absolute path of the manifest, quoted in the refusal
    `short_desc`, `long_desc`    Free text, for humans reading the manifest
    ============================ ===========================================================

    `kinetics_class` is resolved from `kinetics_class_name` against
    :mod:`rmgpy.kinetics` when the manifest is loaded, so a typo raises at database
    load rather than quietly quarantining nothing.
    """

    def __init__(self, family_label, state, kinetics_class_name, reason, path,
                 name='', short_desc='', long_desc=''):
        self.family_label = family_label
        self.state = state
        self.kinetics_class_name = kinetics_class_name
        self.kinetics_class = _resolve_kinetics_class(kinetics_class_name, path)
        self.reason = reason
        self.path = path
        self.name = name
        self.short_desc = short_desc
        self.long_desc = long_desc

    def __repr__(self):
        return '<KineticsQuarantine {0!r} {1} applies to {2}>'.format(
            self.family_label, self.state, self.kinetics_class_name)

    def applies_to(self, kinetics):
        """
        Return ``True`` if `kinetics` falls under this quarantine.

        The test is ``isinstance`` rather than an exact type match, so a future
        subclass of the quarantined model is covered rather than escaping.
        """
        return kinetics is not None and isinstance(kinetics, self.kinetics_class)

    def affected_entries(self, family):
        """
        Enumerate the affected entries of `family` **from the database**, by
        evaluating this manifest's criterion against what is actually loaded.

        Returns a dict with two keys, ``'rules'`` and ``'training'``, each mapping
        to a list of :class:`rmgpy.data.base.Entry`. Nothing here consults a stored
        list of indices or labels, which is what keeps the marker in step with the
        data.
        """
        rules = [entry
                 for entries in getattr(family.rules, 'entries', {}).values()
                 for entry in entries
                 if self.applies_to(entry.data)]

        training = []
        for depository in getattr(family, 'depositories', []):
            if not depository.label.endswith('/training'):
                continue
            training.extend(entry for entry in depository.entries.values()
                            if self.applies_to(entry.data))

        return {'rules': rules, 'training': training}


def _resolve_kinetics_class(name, path):
    """
    Resolve a kinetics class *name* from a manifest against :mod:`rmgpy.kinetics`.

    Raises :class:`DatabaseError` on an unknown name. Failing loudly matters more
    here than anywhere else in this module: a manifest that resolves to nothing
    would report a quarantine that gates nothing, which is worse than no manifest
    at all.
    """
    kinetics_class = getattr(rmgpy.kinetics, name, None)
    if not isinstance(kinetics_class, type):
        raise DatabaseError(
            'Quarantine manifest {0} declares appliesToKineticsClass = {1!r}, which is not a '
            'kinetics class in rmgpy.kinetics. A manifest naming a class that does not exist '
            'would quarantine nothing while claiming to quarantine something.'.format(path, name))
    return kinetics_class


def load_family_quarantine(family_label, family_path):
    """
    Load the quarantine manifest from a kinetics family directory, if it has one.

    Returns a :class:`KineticsQuarantine`, or ``None`` when the family carries no
    manifest -- which is the case for every ordinary family, and is why this costs
    a single :func:`os.path.exists` per family load.

    The manifest is executed the way the rest of the database is, with builtins
    stripped, so it stays a declarative data file rather than a script.
    """
    path = os.path.join(family_path, QUARANTINE_FILENAME)
    if not os.path.exists(path):
        return None

    local_context = {'__builtins__': None}
    global_context = {'__builtins__': None}
    with open(path, 'r') as f:
        content = f.read()
    try:
        exec(content, global_context, local_context)
    except Exception:
        logging.error('Error while reading quarantine manifest %s.', path)
        raise

    missing = [field for field in _REQUIRED_FIELDS if not local_context.get(field)]
    if missing:
        raise DatabaseError(
            'Quarantine manifest {0} is missing required field(s) {1}. A manifest that cannot '
            'say what it applies to and why cannot gate anything.'.format(path, ', '.join(missing)))

    quarantine = KineticsQuarantine(
        family_label=family_label,
        state=local_context['state'],
        kinetics_class_name=local_context['appliesToKineticsClass'],
        reason=local_context['reason'],
        path=path,
        name=local_context.get('name', '') or '',
        short_desc=local_context.get('shortDesc', '') or '',
        long_desc=local_context.get('longDesc', '') or '',
    )
    logging.info('Kinetics family %s is %s (%s): %s',
                 family_label, quarantine.state, quarantine.kinetics_class_name, quarantine.reason)
    return quarantine


def get_quarantine(source):
    """
    Return the :class:`KineticsQuarantine` covering `source`, or ``None``.

    `source` may be a family object or a family/library label. Labels are resolved
    against the loaded kinetics database; an unloaded database, an unknown label, or
    a kinetics library (which cannot carry a family manifest) all yield ``None``.
    """
    if source is None:
        return None
    if hasattr(source, 'quarantine'):
        return source.quarantine

    import rmgpy.data.rmg
    database = getattr(rmgpy.data.rmg, 'database', None)
    if database is None:
        return None
    try:
        family = database.kinetics.families[source]
    except (AttributeError, KeyError, TypeError):
        return None
    return getattr(family, 'quarantine', None)


def describe_provenance(reaction, kinetics=None, source=None, entry=None):
    """
    Describe where a reaction's kinetics came from, for the refusal message.

    `source` and `entry` are the second and third elements of the tuple returned by
    :meth:`KineticsFamily.get_kinetics`, when the caller has them. They are not
    always available: an averaged or generalised rate-rule estimate returns
    ``entry=None``, and callers downstream of estimation (a seed mechanism, a
    reaction library, a resumed model) have neither. The kinetics ``comment``, which
    RMG stamps with the matched node or training reaction, is the fallback, and the
    reaction's own template is the last resort -- so this never returns an empty
    string and the refusal never omits its provenance field.
    """
    parts = []

    if entry is not None:
        depository_label = getattr(source, 'label', None)
        if depository_label:
            parts.append('{0} entry {1} "{2}"'.format(depository_label, entry.index, entry.label))
        else:
            parts.append('rate rule entry {0} "{1}"'.format(entry.index, entry.label))
        if entry.rank is not None:
            parts.append('rank {0}'.format(entry.rank))
    elif isinstance(source, str) and source:
        parts.append('{0} estimate'.format(source))

    template = getattr(reaction, 'template', None)
    if template:
        parts.append('template [{0}]'.format(', '.join(str(label) for label in template)))

    if kinetics is None:
        kinetics = getattr(reaction, 'kinetics', None)
    comment = (getattr(kinetics, 'comment', '') or '').strip()
    if comment:
        parts.append('comment: {0}'.format(' '.join(comment.split())))

    return '; '.join(parts) if parts else 'unrecorded'


def check_quarantine(reaction, stage, family=None, kinetics=None, source=None, entry=None):
    """
    Refuse `reaction` if its kinetics come from quarantined database data.

    Raises :class:`QuarantinedKineticsError` naming the family, the rule or
    training-entry provenance, the generated reaction, the kinetics class, and the
    reason. Returns ``None`` otherwise, including for every reaction of every family
    that carries no manifest -- which is all of ordinary chemistry.

    `stage` names the boundary being crossed, so a reader can tell a refusal at
    kinetics estimation from one at core or edge admission.
    """
    quarantine = get_quarantine(family if family is not None else getattr(reaction, 'family', None))
    if quarantine is None:
        return

    if kinetics is None:
        kinetics = getattr(reaction, 'kinetics', None)
    if not quarantine.applies_to(kinetics):
        return

    try:
        described = str(reaction)
    except Exception:
        described = repr(reaction)

    raise QuarantinedKineticsError(
        '{state}: refusing to admit a quarantined rate at {stage}.\n'
        '  family:             {family}\n'
        '  provenance:         {provenance}\n'
        '  reaction:           {reaction}\n'
        '  kinetics class:     {kinetics_class}\n'
        '  reason:             {reason}\n'
        '  quarantine record:  {path}\n'
        'This reaction is not being dropped, averaged, zeroed, reversed, made irreversible, or '
        'given a substitute rate -- any of those would let the run report success on a mechanism '
        'that is quietly wrong. The family and its data remain in the database untouched. Remove '
        '{path} only when the quarantine is genuinely resolved.'.format(
            state=quarantine.state,
            stage=stage,
            family=quarantine.family_label,
            provenance=describe_provenance(reaction, kinetics=kinetics, source=source, entry=entry),
            reaction=described,
            kinetics_class=type(kinetics).__name__,
            reason=quarantine.reason,
            path=quarantine.path,
        ))
