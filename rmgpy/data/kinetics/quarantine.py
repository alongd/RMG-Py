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

**What is quarantined is the RATE, scoped by the family that authored it.** The defect
is a property of the number: a ``Marcus`` rate with no electrochemical reference returns
a meaningless quantity wherever it is stored, and copying it into a seed mechanism does
not repair it. The family appears in the key for two reasons and neither of them is that
the family is what is wrong -- it is where the manifest is *kept*, next to the data, and
it is what *scopes* the criterion, since ``appliesToKineticsClass = "Marcus"`` on its own
would ban a kinetics class across all of RMG including the legitimate electrochemistry it
was designed for. So the key is (authoring family, kinetics criterion), and "authoring"
is the load-bearing word: it must survive the rate being copied.

That is why the gate resolves the family through :func:`authoring_family` rather than
reading ``reaction.family``. The latter is a slot that two different *kinds* of name
share -- :class:`~rmgpy.data.kinetics.library.LibraryReaction` overwrites it with the
label of the library the rate was loaded from -- so keying on it admits a quarantined
rate copied into a library, and refuses an innocent library whose name collides with a
quarantined family's. Authorship is recovered from the ``family: <label>`` line RMG
writes into an estimated rate's comment and saves into the entry's ``longDesc``. For a
hand-written entry carrying no comment there is nothing to recover, that gap is real,
and :func:`_warn_unattributable` reports it rather than guessing in either direction.

**Every path by which a rate reaches the model, and the gate covering it.** This list is
the answer to "is the gate wired in?", and it is here rather than in a review note so that
a future reviewer can check it against the code. Each entry names the call in
``rmgpy/rmg/model.py``.

======================================= =====================================================
Path                                    Gate
======================================= =====================================================
generated reaction, kinetics estimated  ``apply_kinetics_to_reaction`` -- the primary gate,
                                        upstream of both core and edge, so a refused rate
                                        cannot steer enlargement from the edge either
core admission                          ``add_reaction_to_core`` -- backstop for everything
                                        that never passes through estimation
edge admission                          ``add_reaction_to_edge`` -- same backstop
pressure-dependent network              ``add_reaction_to_unimolecular_networks``. A pdep
                                        path reaction goes *here instead of* core or edge,
                                        and ``generate_kinetics=False`` skips estimation, so
                                        the two backstops above never see it
seed mechanism                          covered: ``add_seed_mechanism_to_core`` reaches
                                        ``add_reaction_to_core`` for every reaction it adds
reaction library                        covered: ``add_reaction_library_to_edge`` reaches
                                        ``add_reaction_to_edge``
library-to-output selection             covered upstream, deliberately not gated again:
                                        ``add_reaction_library_to_output`` only re-selects
                                        reactions already in ``self.edge.reactions``, each of
                                        which passed the edge gate
======================================= =====================================================

What this does **not** cover is unchanged and enumerated per-manifest in ``bypassRoutes``:
an engine without this module, a consumer reading ``.kinetics`` off the database without
model admission at all, and a hand-written library entry carrying no authorship.

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

import ast
import importlib
import inspect
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


def _counts_calls_to(module, symbol_name):
    """
    Return ``(number of calls to `symbol_name` in `module`'s source, source readable)``.

    The bound, stated because the field's name promises more than any static check can
    deliver: this proves a call **exists in the source**, not that it executes on a given
    run or covers every path. What it does close is the gap it was written for -- a module
    that imports the gate and calls it nowhere, which binding alone accepts.
    """
    try:
        tree = ast.parse(inspect.getsource(module))
    except (OSError, TypeError, SyntaxError):
        return 0, False
    calls = 0
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if isinstance(func, ast.Name) and func.id == symbol_name:
            calls += 1
        elif isinstance(func, ast.Attribute) and func.attr == symbol_name:
            calls += 1
    return calls, True


def _check_engine_requirements(path, family_label, manifest):
    """
    Refuse to load a quarantine manifest whose declared engine capability is absent.

    A manifest may declare ``requiresEngineModule`` and ``requiresEngineSymbol``: the module
    and function that make its refusal real. Honouring them here is what turns that
    declaration into a pin rather than a comment -- without it a database can state a
    requirement that nothing reads, and a manifest that is silently inert reads like
    protection while providing none.

    Every field below is either honoured or refused at declaration time. A manifest field
    that is read, found wanting, and then ignored is the defect this function exists to
    remove, so there is no third category:

    * ``requiresEngineModule`` -- must import.
    * ``requiresEngineSymbol`` -- must be present in that module **and callable**. An
      attribute that merely exists proves nothing: ``math.pi`` is not a gate. Declaring a
      symbol without a module is refused rather than ignored, since it asks for a check
      that cannot be performed.
    * ``requiresEngineCallSites`` -- module(s) that must bind the declared symbol to *the
      same object* **and call it**. Binding alone was the first version of this check and
      was not the property the field is named for: deleting every call while keeping the
      import leaves a module that still binds the gate and runs nothing. The call check is
      static (see :func:`_counts_calls_to`), so it proves a call exists in the source, not
      that it executes on a given run.
    * ``requiresEngineCommit`` -- **not declarable**. An installed engine has no reliable
      commit to compare against, so the check would pass on every checkout. Record the
      commit as ``recordedEngineCommit`` instead, which this function reads as provenance
      and never treats as a guarantee.

    One limit remains and cannot be closed from the database side: an engine so old that
    it has no quarantine loader never reaches this function, and never reads the manifest.
    That is the first entry in every manifest's ``bypassRoutes``.
    """
    if manifest.get('requiresEngineCommit') is not None:
        raise DatabaseError(
            'Quarantine manifest {path} (family {family}) declares requiresEngineCommit. That '
            'field is not enforceable: an installed RMG-Py has no reliable commit to compare '
            'against, so the check would pass on every checkout and the manifest would read as '
            'pinned while pinning nothing. Rename it to recordedEngineCommit, which is kept as '
            'provenance, and pin the capability itself with requiresEngineModule, '
            'requiresEngineSymbol and requiresEngineCallSites.'.format(
                path=path, family=family_label))

    module_name = manifest.get('requiresEngineModule')
    symbol_name = manifest.get('requiresEngineSymbol')
    call_sites = manifest.get('requiresEngineCallSites') or ()
    if isinstance(call_sites, str):
        call_sites = (call_sites,)

    if not module_name:
        if symbol_name or call_sites:
            raise DatabaseError(
                'Quarantine manifest {path} (family {family}) declares an engine symbol or call '
                'site without requiresEngineModule, so there is nothing to look them up in and '
                'the requirement cannot be checked. Declare requiresEngineModule too, or remove '
                'the other fields -- an unenforceable requirement must not look like an enforced '
                'one.'.format(path=path, family=family_label))
        return

    try:
        module = importlib.import_module(module_name)
    except ImportError as exc:
        raise DatabaseError(
            'Quarantine manifest {path} requires engine module {module!r}, which this RMG-Py '
            'does not provide ({exc}). The manifest exists because family {family} must not '
            'reach a reaction model on an engine without that capability, so loading is '
            'refused here rather than proceeding unguarded.'.format(
                path=path, module=module_name, exc=exc, family=family_label))

    if not symbol_name:
        if call_sites:
            raise DatabaseError(
                'Quarantine manifest {path} (family {family}) declares requiresEngineCallSites '
                'without requiresEngineSymbol, so there is no symbol whose binding could be '
                'checked at those call sites.'.format(path=path, family=family_label))
        return

    symbol = getattr(module, symbol_name, None)
    if symbol is None:
        raise DatabaseError(
            'Quarantine manifest {path} requires {symbol!r} from engine module {module!r}, '
            'which this RMG-Py provides without it. The capability the manifest depends on '
            'has been renamed or removed, so the refusal it declares for family {family} '
            'would not fire; loading is refused instead.'.format(
                path=path, symbol=symbol_name, module=module_name, family=family_label))

    if not callable(symbol):
        raise DatabaseError(
            'Quarantine manifest {path} requires {symbol!r} from engine module {module!r} as '
            'the gate that refuses family {family}, but this RMG-Py binds that name to '
            '{value!r}, which is not callable. A gate that cannot be called is not a gate, and '
            'accepting any attribute that merely exists would let the pin pass against a '
            'constant.'.format(path=path, symbol=symbol_name, module=module_name,
                               family=family_label, value=symbol))

    for call_site in call_sites:
        try:
            site = importlib.import_module(call_site)
        except ImportError as exc:
            raise DatabaseError(
                'Quarantine manifest {path} requires {symbol!r} to be reached from {site!r}, '
                'which this RMG-Py does not provide ({exc}). The refusal it declares for '
                'family {family} would have no call site.'.format(
                    path=path, symbol=symbol_name, site=call_site, exc=exc,
                    family=family_label))
        if getattr(site, symbol_name, None) is not symbol:
            raise DatabaseError(
                'Quarantine manifest {path} requires {symbol!r} from {module!r} to be wired '
                'into {site!r}, and it is not: {site!r} binds {bound!r} under that name. The '
                'capability exists on this engine but is not reached from the path family '
                '{family} depends on, so the manifest would declare a refusal that never '
                'fires. Existence of the symbol is not evidence that anything calls it, which '
                'is why this is checked separately.'.format(
                    path=path, symbol=symbol_name, module=module_name, site=call_site,
                    bound=getattr(site, symbol_name, None), family=family_label))
        calls, readable = _counts_calls_to(site, symbol_name)
        if not readable:
            raise DatabaseError(
                'Quarantine manifest {path} requires {symbol!r} to be called from {site!r}, '
                'and this engine provides {site!r} without readable source, so the '
                'requirement cannot be checked. Refusing rather than assuming, because an '
                'unverifiable requirement that loads is the failure mode this field group '
                'exists to remove.'.format(path=path, symbol=symbol_name, site=call_site))
        if not calls:
            raise DatabaseError(
                'Quarantine manifest {path} requires {symbol!r} to be reached from {site!r}. '
                '{site!r} imports it and never calls it: {calls:d} call sites in its source. '
                'Deleting every call while keeping the import leaves a module that still '
                'binds the gate and runs nothing, so binding alone is not the property this '
                'field is named for. Family {family} would be admitted unchecked.'.format(
                    path=path, symbol=symbol_name, site=call_site, calls=calls,
                    family=family_label))


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

    _check_engine_requirements(path, family_label, local_context)

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
    recorded_commit = local_context.get('recordedEngineCommit')
    if recorded_commit:
        # Read, so the field is provenance that reaches a reader rather than a line nothing
        # consults. It is deliberately not compared against anything -- see
        # `_check_engine_requirements` for why a commit check could not fail.
        logging.info('  the engine capability it requires was first provided by commit %s '
                     '(provenance only; what is enforced is the capability).', recorded_commit)
    return quarantine


#: Prefix of the line RMG writes into an estimated rate's comment naming the family
#: that authored it. It is written by the kinetics estimator and saved into the entry's
#: ``longDesc`` by the library writer, and is therefore the one piece of authorship that
#: survives a rate being copied into a seed mechanism or a reaction library.
_FAMILY_COMMENT_PREFIX = 'family:'

#: One warning per (library, kinetics class) for rates that match a quarantine criterion
#: but carry no recoverable authorship. Per-reaction logging would emit once per edge
#: reaction per iteration.
_UNATTRIBUTED_WARNED = set()


def _families_from_provenance(*texts):
    """
    Return **every** distinct family label declared by a ``family: <label>`` line.

    Every one, in order of appearance, and not just the first. Provenance is free text:
    taking the first line and stopping means one prepended line shadows the genuine one,
    and a bypass that costs an attacker a single line of comment is not a bound worth
    having. Trusting all of them inverts that -- an added line can only ever *widen* what
    is consulted, never narrow it.

    The trust this places in free text is bounded, and the bound is worth stating. A
    forged or stale line cannot invent a quarantine: the manifest still has to exist, and
    its ``appliesToKineticsClass`` criterion still has to match the rate. What a forged
    line can do is cause a *false refusal* of an independent rate that happens to be of
    the quarantined class. That failure is loud, names the file it came from, and is
    corrected by deleting a line -- where the failure in the other direction is a
    meaningless number in a mechanism that reports success.
    """
    labels = []
    for text in texts:
        for line in (text or '').splitlines():
            stripped = line.strip()
            if stripped.startswith(_FAMILY_COMMENT_PREFIX):
                label = stripped[len(_FAMILY_COMMENT_PREFIX):].strip()
                if label and label not in labels:
                    labels.append(label)
    return labels


def authoring_family(reaction):
    """
    Return the label of the kinetics family that **authored** `reaction`'s rate, or ``None``.

    This is deliberately not ``reaction.family``. That attribute holds two different
    *kinds* of name depending on the wrapper: on a :class:`TemplateReaction` it is the
    family that generated the reaction, but :class:`LibraryReaction` overwrites it with
    the **library** label the reaction was loaded from. Reading it blindly is wrong in
    both directions -- it misses a quarantined rate copied into a library, and it refuses
    an innocent library whose name happens to match a quarantined family's.

    For a library reaction the family is recovered from provenance instead: the
    ``family: <label>`` line RMG writes into an estimated rate's comment and saves into
    the entry's ``longDesc``. That line is genuinely absent from a hand-written entry, in
    which case authorship is unrecoverable and this returns ``None`` -- see
    :func:`check_quarantine` for what is done about that.
    """
    labels = authoring_families(reaction)
    return labels[0] if labels else None


def authoring_families(reaction):
    """
    Return every family label `reaction`'s provenance declares, most likely first.

    :func:`check_quarantine` uses this rather than :func:`authoring_family`, so that a
    second ``family:`` line cannot hide the first. For a template reaction the list is
    the one family that generated it.
    """
    if reaction is None:
        return []
    if getattr(reaction, 'library', None) is None:
        family = getattr(reaction, 'family', None)
        return [family] if family else []
    entry = getattr(reaction, 'entry', None)
    return _families_from_provenance(
        getattr(getattr(reaction, 'kinetics', None), 'comment', ''),
        getattr(entry, 'long_desc', ''),
    )


def iter_quarantines():
    """
    Yield every :class:`KineticsQuarantine` in the loaded kinetics database.

    Empty when no database is loaded, and -- as in all of ordinary chemistry -- when no
    family carries a manifest.
    """
    import rmgpy.data.rmg
    database = getattr(rmgpy.data.rmg, 'database', None)
    families = getattr(getattr(database, 'kinetics', None), 'families', None) or {}
    for family in families.values():
        quarantine = getattr(family, 'quarantine', None)
        if quarantine is not None:
            yield quarantine


def _warn_unattributable(reaction, kinetics):
    """
    Warn when a library rate matches a live quarantine criterion but cannot be attributed.

    Deliberately a warning and not a refusal. The criterion is a kinetics class, which is
    not by itself evidence of quarantined origin -- a ``Marcus`` rate in a genuine
    electrochemistry library is exactly the legitimate case a quarantine must not touch.
    Refusing on the criterion alone would trade a silent admission for a silent
    obstruction. What *is* reportable is the gap itself: a rate that could be a copy of
    quarantined data, in a file that records no authorship to check it against.
    """
    if kinetics is None:
        return
    matches = [q for q in iter_quarantines() if q.applies_to(kinetics)]
    if not matches:
        return
    library = getattr(reaction, 'library', None)
    key = (library, type(kinetics).__name__)
    if key in _UNATTRIBUTED_WARNED:
        return
    _UNATTRIBUTED_WARNED.add(key)
    logging.warning(
        'Reaction library %r supplies %s kinetics, which is the class quarantined for '
        'family %s (%s), and the entry records no authoring family. A rate copied out of '
        'that family would be admitted here unchecked, because nothing in the entry says '
        'where it came from. This is NOT a refusal: the same kinetics class from an '
        'independent source is legitimate and the gate cannot tell the two apart. To '
        'settle it, add a "%s <label>" line to the entry\'s longDesc -- which is what RMG '
        'itself writes for an estimated rate -- or confirm the rate is independent.',
        library, type(kinetics).__name__, matches[0].family_label, matches[0].reason,
        _FAMILY_COMMENT_PREFIX)


#: (database directory, family label) -> KineticsQuarantine or None, for families that
#: are not loaded and had to be answered from disk. A miss costs two `os.path` calls.
_DISK_QUARANTINE_CACHE = {}

#: One warning per (family label, kinetics class) for questions that could not be
#: answered at all.
_UNANSWERED_WARNED = set()


def resolve_quarantine(label):
    """
    Answer "is family `label` quarantined?" as ``(quarantine, answered)``.

    ``answered`` is the whole point of this function. ``get_quarantine`` returns ``None``
    both for "this family carries no manifest" and for "I could not find this family at
    all", and those are opposite situations: the first is a clean bill of health, the
    second is an unanswered question. Collapsing them means a lookup miss silently skips
    every assertion that follows it.

    The unloaded case is not exotic. ``CoreEdgeReactionModel.add_seed_mechanism_to_core``
    deliberately *converts* a seed reaction whose family is unavailable into a library
    reaction rather than loading the family, so "authorship recovered, family not loaded"
    is the normal path for any seed written by a different database than the one running.

    So the question is answered from **disk** when the family is not loaded: a manifest is
    a sidecar file in the family's directory, and reading it needs no family object. What
    is left unanswered after that is the genuine case -- a family this database does not
    contain, whose manifest cannot be consulted because it is not here.
    """
    import rmgpy.data.rmg
    database = getattr(rmgpy.data.rmg, 'database', None)
    families = getattr(getattr(database, 'kinetics', None), 'families', None)
    if families:
        try:
            family = families[label]
        except (KeyError, TypeError):
            pass
        else:
            return getattr(family, 'quarantine', None), True

    from rmgpy import settings
    directory = (settings or {}).get('database.directory')
    if not directory:
        return None, False

    key = (directory, label)
    if key in _DISK_QUARANTINE_CACHE:
        return _DISK_QUARANTINE_CACHE[key]

    family_path = os.path.join(directory, 'kinetics', 'families', str(label))
    if not os.path.isdir(family_path):
        answer = (None, False)
    elif not os.path.exists(os.path.join(family_path, QUARANTINE_FILENAME)):
        # The family is in this database and carries no manifest. That is an answer.
        answer = (None, True)
    else:
        answer = (load_family_quarantine(label, family_path), True)
    _DISK_QUARANTINE_CACHE[key] = answer
    return answer


def _warn_unanswered(reaction, kinetics, labels):
    """
    Report a rate whose authorship is known and whose manifest could not be consulted.

    Deliberately not a refusal, and the reasoning is the opposite of
    :func:`_warn_unattributable`'s. There the label was missing; here the family is. A
    seed mechanism that names families this database does not contain is ordinary -- the
    reaction model converts such reactions on purpose and merely logs it -- so refusing
    would stop runs that have nothing to do with any quarantine, and would break the
    promise that a database with no manifest anywhere behaves exactly as it did before.

    What must not happen is that it passes in silence, which is what a bare ``None``
    achieved. The residual risk is stated rather than closed: if the missing family is
    quarantined *in its own database*, this run cannot know it.
    """
    if kinetics is None:
        return
    if not any(True for _ in iter_quarantines()):
        # No quarantine anywhere in this database: nothing is being skipped, and warning
        # here would put a line into every ordinary run that loads a foreign seed.
        return
    key = (tuple(labels), type(kinetics).__name__)
    if key in _UNANSWERED_WARNED:
        return
    _UNANSWERED_WARNED.add(key)
    logging.warning(
        'Cannot tell whether %s kinetics authored by family %s are quarantined: this '
        'database has no such family, so its manifest -- if it has one -- cannot be '
        'consulted. The rate IS being admitted. This is not the same as a rate with no '
        'recorded authorship: the authorship is here and the answer is not. If that '
        'family is quarantined in the database it came from, this run cannot know it. '
        'Load that family, or supply the rate from a source this database can check.',
        type(kinetics).__name__, ', '.join(str(label) for label in labels))


def get_quarantine(source):
    """
    Return the :class:`KineticsQuarantine` covering `source`, or ``None``.

    `source` may be a family object or a family label. Labels are resolved against the
    loaded kinetics database; an unloaded database and an unknown label both yield
    ``None``. Callers must not pass a *library* label here -- see
    :func:`authoring_family` for why the two kinds of name are not interchangeable.

    **Prefer :func:`resolve_quarantine` in a gate.** This function cannot distinguish
    "no manifest" from "no such family", and a gate that cannot tell those apart admits
    on a lookup miss.
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

    The family is resolved by :func:`authoring_family`, which asks who *authored* the
    rate rather than reading ``reaction.family`` -- a slot whose meaning changes with the
    wrapper. When a caller passes `family` explicitly it is used as given, because the
    caller is then holding the family object itself.
    """
    if kinetics is None:
        kinetics = getattr(reaction, 'kinetics', None)

    if family is not None:
        # The caller is holding the family itself; nothing to resolve or doubt.
        candidates = [(get_quarantine(family), True, family)]
    else:
        labels = authoring_families(reaction)
        if not labels:
            if getattr(reaction, 'library', None) is not None:
                _warn_unattributable(reaction, kinetics)
            return
        candidates = [resolve_quarantine(label) + (label,) for label in labels]

    unanswered = [label for quarantine, answered, label in candidates if not answered]
    for quarantine, answered, _label in candidates:
        if answered and quarantine is not None and quarantine.applies_to(kinetics):
            break
    else:
        # Nothing matched. Say so where the question could not be asked, rather than
        # letting a lookup miss look like a clean bill of health.
        if unanswered:
            _warn_unanswered(reaction, kinetics, unanswered)
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
