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
import stat

import rmgpy.kinetics
from rmgpy.exceptions import DatabaseError, QuarantinedKineticsError

#: Filename of the sidecar manifest inside a kinetics family directory.
QUARANTINE_FILENAME = 'quarantine.py'

#: Manifest fields without which a quarantine cannot say what it covers or why.
#: ``name``, ``shortDesc`` and ``longDesc`` are optional and free text; they follow the
#: rest of the database in being camelCase in the data file and snake_case on the object.
#: Distinguishes a manifest field DECLARED empty from one not declared at all.
_ABSENT = object()

_REQUIRED_FIELDS = ('state', 'appliesToKineticsClass', 'reason')


#: The admission paths this engine defines, per module: every function through which a
#: reaction can enter a model without passing one of the others. A manifest that names a
#: module here is asking for *those* sites to be gated, not for a call to appear somewhere
#: in the file. Keeping this list beside the checker is deliberate -- the manifest cannot
#: know an engine's internal structure, and an engine that grows a fifth admission path
#: should fail its own pin until the path is gated.
_GATE_CALL_SITES = {
    'rmgpy.rmg.model': (
        'apply_kinetics_to_reaction',            # estimation
        'add_reaction_to_core',                  # seeds, libraries, API callers
        'add_reaction_to_edge',                  # the edge half of the same
        'add_reaction_to_unimolecular_networks',  # pdep, which reaches neither of those
    ),
}


def _is_statically_dead(node):
    """Whether `node` is an `if`/`while` whose test can never be true."""
    test = getattr(node, 'test', None)
    if test is None:
        return False
    if isinstance(test, ast.Constant):
        return not test.value
    return False


def _gate_calls_in_source(module, symbol_name, module_name):
    """
    Return ``({enclosing function -> call count}, source readable)``: a SYNTACTIC census.

    **This is a static syntactic presence check and nothing more.** The name says exactly
    that, after four consecutive review rounds each found the previous name and docstring
    claiming more than the code delivers -- the last of them for the phrase this paragraph
    replaces, which said the call was in a "reachable branch". It is not a reachability
    analysis. What it reports is: for each function named in :data:`_GATE_CALL_SITES`, how
    many ``ast.Call`` nodes inside that function's body call `symbol_name`, either by name
    or through an alias bound to `module_name`.

    Three shapes are excluded, because each let a module satisfy the pin while gating
    nothing and each is cheap to spot in the tree:

    * calls inside an ``if``/``while`` whose test is a falsy literal;
    * calls at module scope rather than inside a function, which run once at import and
      gate no reaction;
    * ``something.check_quarantine(...)`` where ``something`` is not the declared engine
      module -- a same-named method on an unrelated object.

    What it does NOT verify, enumerated rather than gestured at, because an incomplete
    list is how the overstatement kept coming back:

    * that the call **executes** on any given run;
    * that it **dominates** the admission it is supposed to guard -- a call placed after
      the reaction has already been appended satisfies this check;
    * that its **arguments** are the reaction being admitted;
    * that the exception it raises **propagates** rather than being swallowed;
    * short-circuited calls such as ``False and check_quarantine(r)``, which this accepts;
    * a **locally shadowed** ``check_quarantine`` bound to something else inside the
      function;
    * calls hidden in lambdas, comprehensions or nested definitions that are never called.

    Deliberately not strengthened again. Each narrowing removes one more example while
    leaving the class untouched, and makes the check *look* stronger without changing what
    it can promise; a reader who trusts a syntactic check to prove execution is misled by
    a better one exactly as much as by this one. Execution is provable, but by coverage
    instrumentation on a real run, in CI -- not from a manifest.

    What it does close, and the reason it is worth having: **coverage of the named sites.**
    With :data:`_GATE_CALL_SITES`, deleting the gate from one admission path out of four is
    refused, where a module-level call count accepted it, and an engine that grows a fifth
    admission path fails its own pin until that path is gated.
    """
    try:
        tree = ast.parse(inspect.getsource(module))
    except (OSError, TypeError, SyntaxError):
        return {}, False

    # Names bound to the engine module itself, so `q.check_quarantine(...)` counts while
    # `somebody_else.check_quarantine(...)` does not.
    module_aliases = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == module_name:
                    module_aliases.add(alias.asname or alias.name.split('.')[0])

    found = {}

    def visit(node, function):
        for child in ast.iter_child_nodes(node):
            if _is_statically_dead(child):
                # Walk the else branch only: the body is unreachable.
                for orelse in getattr(child, 'orelse', []):
                    visit(orelse, function)
                continue
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
                visit(child, child.name)
                continue
            if isinstance(child, ast.Call) and function is not None:
                func = child.func
                if isinstance(func, ast.Name) and func.id == symbol_name:
                    found[function] = found.get(function, 0) + 1
                elif (isinstance(func, ast.Attribute) and func.attr == symbol_name
                      and isinstance(func.value, ast.Name)
                      and func.value.id in module_aliases):
                    found[function] = found.get(function, 0) + 1
            visit(child, function)

    visit(tree, None)
    return found, True


#: One message per manifest still using the pre-round-99 spelling of the field.
_LEGACY_CALL_SITES_WARNED = set()


def _warn_legacy_call_sites_field(path):
    """
    Report a manifest declaring ``requiresEngineCallSites``, once per manifest.

    The field was renamed because its name claimed more than the check delivers -- it is a
    static syntactic presence check, not a verification of call sites. The old spelling is
    still read rather than ignored: manifests declaring it ship in database repositories
    that this engine change cannot edit, and quietly dropping the field would disarm a pin
    in a repository nobody is looking at. A warning, not a refusal, for the same reason.
    """
    if path in _LEGACY_CALL_SITES_WARNED:
        return
    _LEGACY_CALL_SITES_WARNED.add(path)
    logging.warning(
        'Quarantine manifest %s declares requiresEngineCallSites. That field is now '
        'requiresEngineCallsInSource, which says what is actually checked: that the '
        'named modules bind the gate and that their SOURCE TEXT contains a call to it. '
        'It is a static syntactic check and never evidence that the call executes. The '
        'old spelling is still honoured; rename it in the database when convenient.', path)


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
    * ``requiresEngineCallsInSource`` -- module(s) that must bind the declared symbol to
      *the same object* and whose **source text contains a call to it** inside each of the
      functions named in :data:`_GATE_CALL_SITES`. That is the whole of it: a **static
      syntactic presence check**. It does not show the call executes, dominates the
      admission it guards, receives the right arguments, or propagates what it raises --
      :func:`_gate_calls_in_source` enumerates the shapes that satisfy it while gating
      nothing, and five review rounds have each found one more. Binding alone is still
      weaker and still refused: deleting every call while keeping the import leaves a
      module that binds the gate and runs nothing.

      Spelled ``requiresEngineCallSites`` until round 99. That name claimed the field
      verified call *sites*, which is a property of a running program, and the ledger
      inherited the claim. The old spelling is still honoured, with a warning, because
      manifests that declare it ship in a database repository this change may not edit;
      declaring both is refused rather than silently resolved.
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
            'requiresEngineSymbol and requiresEngineCallsInSource.'.format(
                path=path, family=family_label))

    module_name = manifest.get('requiresEngineModule')
    symbol_name = manifest.get('requiresEngineSymbol')
    # Presence, not truthiness. `manifest.get(...) or ()` collapses "declared empty" into
    # "not declared", and both promises made here are about DECLARING a field: an empty
    # new field beside a populated old one slipped through the ambiguity refusal, which is
    # precisely the manifest someone writes halfway through a rename.
    declared = manifest.get('requiresEngineCallsInSource', _ABSENT)
    legacy = manifest.get('requiresEngineCallSites', _ABSENT)
    if declared is not _ABSENT and legacy is not _ABSENT:
        raise DatabaseError(
            'Quarantine manifest {path} (family {family}) declares both '
            'requiresEngineCallsInSource and its old spelling requiresEngineCallSites. '
            'Resolving that quietly would mean guessing which one the author meant to be '
            'read; declare one.'.format(path=path, family=family_label))
    if legacy is not _ABSENT:
        _warn_legacy_call_sites_field(path)
        declared = legacy

    if declared is _ABSENT:
        call_sites = ()
    else:
        call_sites = (declared,) if isinstance(declared, str) else tuple(declared)
        if not call_sites:
            raise DatabaseError(
                'Quarantine manifest {path} (family {family}) declares a call-in-source '
                'requirement that names no module, so it pins nothing while reading as a '
                'pin. Name the modules that must call the gate, or remove the field.'.format(
                    path=path, family=family_label))

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
                'Quarantine manifest {path} (family {family}) declares '
                'requiresEngineCallsInSource without requiresEngineSymbol, so there is no '
                'symbol whose calls could be looked for in those modules.'.format(
                    path=path, family=family_label))
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
        gated, readable = _gate_calls_in_source(site, symbol_name, module_name)
        calls = sum(gated.values())
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

        required = _GATE_CALL_SITES.get(call_site)
        if required:
            missing = [name for name in required if name not in gated]
            if missing:
                raise DatabaseError(
                    'Quarantine manifest {path} requires {symbol!r} to be reached from '
                    '{site!r}, and {site!r} gates only {present}. This engine admits a '
                    'reaction through {total:d} paths in that module and {count:d} of them '
                    'call the gate: {missing} do not. A rate entering by one of those is '
                    'admitted without the manifest for family {family} ever being consulted, '
                    'so a module-level call count is not the property this field is named '
                    'for.'.format(
                        path=path, symbol=symbol_name, site=call_site,
                        present=', '.join(sorted(gated)) or 'nothing',
                        total=len(required), count=len(required) - len(missing),
                        missing=', '.join(missing), family=family_label))


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


#: One message per family directory whose manifest was refused as a file.
_UNSAFE_MANIFESTS_WARNED = set()


def _warn_unsafe_manifest(family_path, reason):
    """Report a manifest refused as a file, once per family directory. Loud on purpose."""
    if family_path in _UNSAFE_MANIFESTS_WARNED:
        return
    _UNSAFE_MANIFESTS_WARNED.add(family_path)
    logging.error(
        'Refusing to read the quarantine manifest in %s: %s. A manifest is a regular file '
        'inside its own family directory, and it is EXECUTED -- so one that is a link to '
        'somewhere else, or is not a file at all, is refused rather than followed, and '
        'one that cannot be examined at all is refused rather than assumed away. The '
        'family is treated as unanswered, not as unquarantined.', family_path, reason)


def _anchor_and_components(family_path):
    """
    Split `family_path` into ``(anchor, components)`` for a link-free descent, or
    ``(None, None)`` when it cannot be split into one.

    The anchor is the one directory that is opened **following** links, because a path
    cannot be descended without something to descend from and that first open always
    resolves whatever it is given. Every component below it is opened with ``O_NOFOLLOW``,
    so the rule is exactly one sentence: nothing below the anchor may be a link.

    The anchor is ``settings['database.directory']`` whenever `family_path` lies below it
    lexically -- which puts the followed open on *local configuration* and pins everything
    a shipped ``family:`` line can steer, including ``kinetics`` and ``families``. Round 99
    anchored on the family's own parent instead, so ``kinetics/families`` was followed; a
    link there escapes the database entirely, and ``_family_directory``'s containment test
    cannot see it, because with that link in place the root and the family resolve into
    the same foreign tree and the prefix comparison passes.

    Falls back to the family's parent when the path is not below the configured database
    -- a family loaded from somewhere else, or a caller passing its own directory. That
    caller is vouching for the path it passed; what is still pinned is the family
    directory and the manifest.

    Refuses outright, rather than approximating, when the split would yield an empty
    component, ``.`` or ``..``. ``os.path.split('')`` gives ``('', '')`` and
    ``os.path.join('', 'quarantine.py')`` is a *relative* name, so the old fallback read
    ``quarantine.py`` out of the process's working directory.
    """
    from rmgpy import settings

    if not family_path:
        # `os.path.relpath('')` raises, and a gate that raises out of the path it is
        # refusing has failed in the direction this whole module exists to avoid.
        return None, None

    root = (settings or {}).get('database.directory')
    components = None
    if root:
        try:
            relative = os.path.relpath(family_path, root)
        except (ValueError, TypeError):
            relative = os.pardir
        if not os.path.isabs(relative) and relative.split(os.sep)[0] != os.pardir:
            anchor, components = root, relative.split(os.sep)

    if components is None:
        anchor, name = os.path.split(family_path.rstrip(os.sep) or family_path)
        components = [name]

    if not anchor or any(part in ('', os.curdir, os.pardir) for part in components):
        return None, None
    return anchor, components


def _read_manifest(family_path):
    """
    Return the manifest's text from inside `family_path`, or ``None``.

    The file at the end of this path is **executed**, and that is what sets the standard
    for how it may be addressed. Round 92 validated the family *directory* and then opened
    ``<directory>/quarantine.py`` by name. Round 95 closed the manifest half of that by
    reading the file through a descriptor on its directory -- and left the directory itself
    named twice: :func:`_family_directory` resolves it with ``realpath`` and returns the
    *unresolved* join, and this function then opened that string with ``O_DIRECTORY`` and
    no ``O_NOFOLLOW``. Check and use were two independent resolutions of one name, so a
    directory symlink swapped into the window between them was followed all the way to
    ``exec()``. Round 95's claim here that "the descriptor pins the directory that was
    validated" was simply false: a descriptor pins whatever the name resolved to when it
    was opened, which is not the same thing as what was checked.

    What closes it is structure, not a third check. The label is a single path component
    by the time it arrives -- :func:`_family_directory` refuses separators -- so the whole
    descent is two steps: open the parent, open the label through it with ``O_NOFOLLOW``,
    open the manifest through *that* with ``O_NOFOLLOW``. Neither component can be a link
    and neither is re-derived from a string after being examined, so "the executed file is
    a direct child of a direct child of the families root" holds by construction and has
    no window in it. The parent is deliberately still followed: it comes from
    ``settings['database.directory']``, which is local configuration, and symlinking a
    database checkout is ordinary. The label is the part a shipped ``family:`` line
    controls, and the label is the part that is pinned.

    The manifest must also be a regular file: a FIFO would otherwise block a run forever
    on the read, and a directory would raise from inside a gate whose entire purpose is to
    keep running and refuse.

    ``None`` covers three cases the caller separates by other means: no manifest (every
    ordinary family), a manifest refused as a file, and one that vanished mid-read.
    """
    if os.open not in os.supports_dir_fd:
        # Without directory descriptors no component can be opened relative to one, so
        # there is no way to descend without following links and no containment to offer.
        # Round 99 fell back to opening the joined pathname and answering anyway, with a
        # comment admitting the hole. A documented hole is still a hole: this module's
        # whole purpose is to refuse rather than to answer unsafely, and it refuses here.
        _warn_unsafe_manifest(
            family_path,
            'this platform provides no directory descriptors, so no component of the '
            'path can be opened without following links')
        return None, None

    anchor, components = _anchor_and_components(family_path)
    if anchor is None:
        _warn_unsafe_manifest(
            family_path,
            'it does not name a family directory below a directory this process may '
            'anchor on')
        return None, None

    nofollow = getattr(os, 'O_NOFOLLOW', 0)
    directory = getattr(os, 'O_DIRECTORY', 0)
    descriptors = []
    fd = None
    try:
        try:
            descriptors.append(os.open(anchor, os.O_RDONLY | directory))
            for component in components:
                descriptors.append(os.open(component, os.O_RDONLY | directory | nofollow,
                                           dir_fd=descriptors[-1]))
            fd = os.open(QUARANTINE_FILENAME, os.O_RDONLY | nofollow,
                         dir_fd=descriptors[-1])
        except (FileNotFoundError, NotADirectoryError):
            return None, None
        except OSError as error:
            # ELOOP lands here for every link below the anchor -- an intermediate
            # directory, the family directory, or the manifest itself. All are refused
            # rather than resolved, and the caller turns the refusal into an unanswered
            # question, never a clean bill.
            _warn_unsafe_manifest(family_path, error)
            return None, None

        info = os.fstat(fd)
        if not stat.S_ISREG(info.st_mode):
            _warn_unsafe_manifest(family_path, 'it is not a regular file')
            return None, None
        # The identity comes from THIS descriptor, not from a second look at the name.
        # `_manifest_signature` stats the path independently, and a rename between that
        # stat and this open would otherwise cache the content of B under the identity of
        # A -- after which restoring A produces a cache HIT returning B's quarantine, and
        # A's criterion is silently bypassed. `st_ctime_ns` cannot help: it would be A's.
        identity = (info.st_mtime_ns, info.st_ctime_ns, info.st_size, info.st_ino)
        with os.fdopen(fd, 'r') as handle:
            fd = None                    # fdopen owns it now
            return handle.read(), identity
    except OSError:
        return None, None
    finally:
        for descriptor in [fd] + descriptors:
            if descriptor is not None:
                try:
                    os.close(descriptor)
                except OSError:
                    pass


def load_family_quarantine(family_label, family_path):
    """
    Load the quarantine manifest from a kinetics family directory, if it has one.

    Returns a :class:`KineticsQuarantine`, or ``None`` when the family carries no
    manifest -- which is the case for every ordinary family, and is why this costs
    a single failed :func:`os.open` per family load.

    The manifest is executed the way the rest of the database is, with builtins
    stripped, so it stays a declarative data file rather than a script.
    """
    return _load_with_identity(family_label, family_path)[0]


def _load_with_identity(family_label, family_path):
    """
    :func:`load_family_quarantine`, also returning the identity of the file it parsed.

    The identity is taken with :func:`os.fstat` on the descriptor the bytes were read
    from, so it names the file that was actually executed rather than whatever the path
    resolved to when somebody else stat'd it. :func:`resolve_quarantine` caches under
    this, which is what makes the cached answer and the cached key describe one object.
    """
    path = os.path.join(family_path, QUARANTINE_FILENAME)
    content, identity = _read_manifest(family_path)
    if content is None:
        # No manifest, or one this database may not execute, or one that vanished between
        # being seen and being read. All three are ``None`` here and the caller decides
        # what that means: `resolve_quarantine` knows which of the three it saw, from the
        # signature it took, and turns the last two into an unanswered question rather
        # than a clean bill of health. Raising instead would kill a run over a concurrent
        # database edit.
        return None, None

    local_context = {'__builtins__': None}
    global_context = {'__builtins__': None}
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
    return quarantine, identity


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
    second ``family:`` line cannot hide the first.

    That guarantee used to hold for a library reaction only. A template reaction took a
    short path -- return the single ``reaction.family`` slot -- and that slot is written
    once per ``family:`` line by :meth:`KineticsLibrary.get_library_reactions`, so only
    the last one survived. An entry naming a quarantined family and then an ordinary
    loaded one was therefore admitted: the last label named a family that *is* loaded, so
    ``CoreEdgeReactionModel`` never converted the reaction to a library reaction, and the
    short path never reached the entry round 92 had just attached. The carrier was there;
    the reader did not read it.

    So provenance is consulted for every shape, and the slot is *added to* what provenance
    gives rather than replacing it. The result is a superset of what the slot alone said,
    which is the direction that matters: an extra label can cause a false refusal, which
    is loud and is corrected by deleting a line, where a missing one admits a rate the
    database has refused.
    """
    if reaction is None:
        return []
    entry = getattr(reaction, 'entry', None)
    labels = _families_from_provenance(
        getattr(getattr(reaction, 'kinetics', None), 'comment', ''),
        getattr(entry, 'long_desc', ''),
    )
    if getattr(reaction, 'library', None) is None:
        # A template reaction's `family` slot really is authorship. A library reaction's
        # is not -- the same slot holds the LIBRARY label there (`LibraryReaction.__init__`),
        # and `electron_placement.py` reads it that way deliberately -- so it is only safe
        # to add for this shape.
        slot = getattr(reaction, 'family', None)
        if slot and slot not in labels:
            labels.append(slot)
    return labels


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


#: families directory -> (signature, whether any manifest was found there).
_DISK_ANY_QUARANTINE_CACHE = {}


def _database_has_any_quarantine():
    """
    Whether this database quarantines anything at all, loaded **or on disk**.

    This is the suppression test for :func:`_warn_unanswered`, and it used to iterate
    loaded families only. That made the warning disappear exactly when it was needed:
    a run early enough that nothing has been loaded yet, in a database that does carry
    manifests, admitted a rate whose authorship it could not resolve and said nothing.
    The policy that admission is allowed *because* it is loud cannot rest on a test that
    goes quiet first.

    Bound, stated because it is real: the disk half is keyed on the families directory's
    own signature, so a manifest dropped into an already-existing family directory during
    a run is not noticed here until something re-reads that directory. It can only make
    this warning appear or not appear -- never a refusal, which goes through
    :func:`resolve_quarantine` and re-checks the manifest's own signature every time.
    """
    if any(True for _ in iter_quarantines()):
        return True

    from rmgpy import settings
    directory = (settings or {}).get('database.directory')
    if not directory:
        return False
    families_root = os.path.join(directory, 'kinetics', 'families')
    try:
        info = os.stat(families_root)
    except OSError:
        return False
    signature = (info.st_mtime_ns, info.st_ino)

    cached = _DISK_ANY_QUARANTINE_CACHE.get(families_root)
    if cached is not None and cached[0] == signature:
        return cached[1]

    found = False
    try:
        names = os.listdir(families_root)
    except OSError:
        names = []
    for name in names:
        if os.path.lexists(os.path.join(families_root, name, QUARANTINE_FILENAME)):
            found = True
            break
    if found:
        # Only the positive is cached, and the asymmetry is the repair. The key is the
        # families root's own signature, and writing a manifest *inside* an existing
        # family directory does not move it -- so a cached negative outlived the fact it
        # recorded, and `_warn_unanswered` stayed silent in a database that had started
        # quarantining something. A cached positive cannot fail that way: the worst it can
        # do is keep a warning switched on after the last manifest is deleted, and this
        # suppression only ever decides whether a warning is printed, never whether a rate
        # is refused. The re-scan a missing negative costs is one `listdir` plus one
        # `lexists` per family, and it happens once per distinct unanswered question --
        # `_warn_unanswered` de-duplicates before it asks.
        _DISK_ANY_QUARANTINE_CACHE[families_root] = (signature, True)
    return found


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


#: (database directory, family label) -> (signature, answer), for families that are not
#: loaded and had to be answered from disk. The signature is what makes this a cache and
#: not a snapshot: it is re-taken on every lookup -- two `os.stat` calls -- and the stored
#: answer is used only while it still matches. Caching the answer alone, as the first
#: version did, meant a manifest added, removed or edited after the first lookup was never
#: seen again for the life of the process, in either direction.
_DISK_QUARANTINE_CACHE = {}


def _family_directory(families_root, label):
    """
    Resolve `label` to a directory inside `families_root`, or ``None``.

    A family label is not a path. It arrives from a ``family:`` line in an entry's
    ``longDesc`` -- user-authored text, in a file that may have been written by another
    database, another group, or an attacker -- and the file at the end of the path this
    function builds is then **executed** by :func:`load_family_quarantine`. Passing that
    string to ``os.path.join`` unexamined meant an absolute label discarded the database
    prefix entirely, and ``..`` segments walked out of the database.

    The checks here are **diagnostic**, and it is worth being exact about that, because
    round 95 mistook them for the closure and left a race behind. Every one of them
    examines a *name*; the file is opened from that name afterwards, and between the two
    the filesystem can change. What actually keeps the executed file inside the database
    is :func:`_read_manifest` opening the label as one ``O_NOFOLLOW`` component below the
    families root, which cannot be raced because it never resolves the name twice. These
    checks exist to refuse early and to say *why* in a message a database author can act
    on -- "it contains a path separator" is worth a great deal more than ``ELOOP``.

    Three of them. The first is syntactic: the label must be a single path component, so
    no separator, no ``..``, no drive letter, nothing absolute, and no control character.
    The second refuses a family directory that is a symbolic link, of any kind, pointing
    anywhere -- a family directory is a real directory, a direct child of the families
    root, which is a rule with no window in it. The third is positional and catches what
    syntax cannot, by requiring the resolved real path to sit under the resolved real
    families root; it is kept because it names the escape precisely when an intermediate
    component is the link.

    Returns the joined path (which need not exist; the caller distinguishes "no such
    family" from "no manifest") or ``None`` when the label is not a name this database
    could have given a family.
    """
    text = str(label)
    if not text or text in ('.', '..'):
        _warn_unsafe_label(label, 'it is empty or names a directory relative to itself')
        return None
    control = next((character for character in text if character < ' '), None)
    if control is not None:
        # A NUL reached `os.stat` and came back out as an uncaught ValueError, killing the
        # run from inside a gate whose whole purpose is to keep running and refuse. A
        # newline is worse than noisy: a directory really can be named with one, so the
        # label selected a family and its manifest was executed. No family this database
        # could have written has a control character in its name.
        _warn_unsafe_label(label, 'it contains the control character {0!r}'.format(control))
        return None
    if os.path.isabs(text) or os.path.splitdrive(text)[0]:
        _warn_unsafe_label(label, 'it is an absolute path')
        return None
    if os.path.basename(text) != text:
        _warn_unsafe_label(label, 'it contains a path separator')
        return None

    family_path = os.path.join(families_root, text)
    if os.path.islink(family_path):
        _warn_unsafe_label(label, 'the family directory is a symbolic link')
        return None
    try:
        resolved_root = os.path.realpath(families_root)
        resolved_family = os.path.realpath(family_path)
    except OSError:
        return None
    if resolved_family != resolved_root and not resolved_family.startswith(
            resolved_root + os.sep):
        _warn_unsafe_label(label, 'it resolves to {0}, outside the database'.format(
            resolved_family))
        return None
    return family_path


def _warn_unsafe_label(label, reason):
    """
    Report a family label that was refused as a path, once per label.

    Loud on purpose. An ordinary database never produces one of these, so either a
    library was written by something that does not understand the format, or the
    ``family:`` line is hostile -- and the file it points at would have been executed.
    """
    if label in _UNSAFE_LABELS_WARNED:
        return
    _UNSAFE_LABELS_WARNED.add(label)
    logging.error(
        'Refusing to look for a quarantine manifest under the family label %r: %s. A '
        'family label names a directory inside the database\'s kinetics/families, and '
        'the manifest found there is executed. Rates attributed to this label are '
        'treated as having an unresolvable authorship.', label, reason)


#: One message per label refused by :func:`_family_directory`.
_UNSAFE_LABELS_WARNED = set()


def _manifest_signature(family_path):
    """
    ``(family directory exists, manifest identity or None, what it is)`` for `family_path`.

    The identity is ``(mtime_ns, ctime_ns, size, inode)``. ``ctime_ns`` is there because
    the other three are all restorable: an in-place edit of the same length, followed by
    ``os.utime`` putting the modification time back, leaves ``(mtime_ns, size, inode)``
    byte-identical while the file says something different -- and the stale answer then
    survives for the life of the run. No userspace call can set the inode change time, so
    it moves on any write. It costs nothing; the ``lstat`` is already being taken.

    Hashing the content would be stronger and is deliberately not done: it turns every
    lookup of every family into a file read, to cover a case ``ctime`` already covers
    anywhere ``ctime`` is real. Where it is not -- a filesystem that reports it as a copy
    of ``mtime`` -- this is no weaker than it was.

    ``kind`` is one of four, and the last is the point of the function:

    * ``'absent'`` -- there is genuinely no manifest here. Every ordinary family.
    * ``'regular'`` -- a manifest that may be read.
    * ``'other'`` -- something is there and it is not a file this database may execute.
    * ``'unreadable'`` -- a manifest **may or may not** be there and this process could not
      find out: a permission error, an I/O error, a descriptor exhaustion. Round 99's
      defect was mapping this onto ``'absent'``, which the caller turns into "this family
      carries no manifest" -- a clean bill of health produced by the failure of the check
      itself, and then cached. It is round 89's collapse of "no" into "I could not look",
      one layer down, and the two are kept apart here for the same reason.

    Taken with :func:`os.lstat` rather than :func:`os.stat`, so that a manifest which is a
    *link* is seen as a link and reported as ``'other'`` instead of as an ordinary file.
    Retargeting the link moves the link's own timestamps, so the signature still
    invalidates.
    """
    try:
        info = os.lstat(os.path.join(family_path, QUARANTINE_FILENAME))
    except (FileNotFoundError, NotADirectoryError):
        return os.path.isdir(family_path), None, 'absent'
    except ValueError:
        # An embedded NUL cannot name a file on any filesystem, so there is nothing here
        # to be unsure about. `_family_directory` refuses these by name; this is a
        # backstop for the paths that do not come through it.
        return os.path.isdir(family_path), None, 'absent'
    except OSError:
        return os.path.isdir(family_path), None, 'unreadable'
    identity = (info.st_mtime_ns, info.st_ctime_ns, info.st_size, info.st_ino)
    kind = 'regular' if stat.S_ISREG(info.st_mode) else 'other'
    return os.path.isdir(family_path), identity, kind

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

    A *loaded* family goes through the same disk read, and used to not. Returning
    ``family.quarantine`` straight from the object meant a manifest added, edited or
    removed after the database was loaded was invisible for the life of the process --
    which is precisely the staleness round 92 repaired on the unloaded half and left
    standing on this one. One path, signature-validated, for both.
    """
    import rmgpy.data.rmg
    database = getattr(rmgpy.data.rmg, 'database', None)
    families = getattr(getattr(database, 'kinetics', None), 'families', None)
    loaded = None
    if families:
        try:
            loaded = families[label]
        except (KeyError, TypeError):
            loaded = None

    from rmgpy import settings
    directory = (settings or {}).get('database.directory')
    if not directory:
        # Nothing to re-read against. A loaded family's own attribute is then the only
        # answer there is, and it is still an answer.
        if loaded is not None:
            return getattr(loaded, 'quarantine', None), True
        return None, False

    families_root = os.path.join(directory, 'kinetics', 'families')
    family_path = getattr(loaded, 'quarantine_path', None) or _family_directory(
        families_root, label)
    if family_path is None:
        # The label did not name a family in this database's families directory. It is
        # not an answer, and it is emphatically not licence to read the file it points
        # at: see _family_directory.
        return None, False

    key = (directory, label)
    signature = _manifest_signature(family_path)
    cached = _DISK_QUARANTINE_CACHE.get(key)
    if cached is not None and cached[0] == signature:
        return cached[1]

    directory_exists, manifest_stat, kind = signature
    if kind == 'unreadable':
        # FIRST, before existence. `directory_exists` is `os.path.isdir`, which returns
        # False when it cannot look -- so a permission error on `kinetics/families` makes
        # the family read as "not where this database would put it", and a LOADED family
        # then answered from its own attribute with answered=True. For a family whose
        # attribute is None that is a clean bill of health produced by the failure of the
        # check itself: round 99's defect, surviving on the branch round 99 did not test.
        #
        # There may or may not be a manifest here; this process could not find out. It is
        # an unanswered question -- and returned UNCACHED, because the condition that
        # produced it (a permission, an I/O error, a descriptor exhaustion) is transient
        # by nature and the next call may well succeed. Caching it would make a momentary
        # failure permanent for the life of the run.
        _warn_unsafe_manifest(family_path, 'it could not be examined')
        return None, False
    if not directory_exists:
        if loaded is not None:
            # A loaded family whose directory is not where this database would put it --
            # a family loaded from elsewhere, or a test's stub -- still knows what it read
            # at load time. Returned UNCACHED, deliberately: this cache is keyed on a
            # directory and a label, and an answer that came from a family *object* is not
            # a property of that key. Caching it made the next lookup of the same label,
            # with no database loaded at all, inherit a quarantine from an object that was
            # no longer there.
            return getattr(loaded, 'quarantine', None), True
        answer = (None, False)
    elif manifest_stat is None:
        # The family is in this database and carries no manifest. That is an answer.
        answer = (None, True)
    elif kind != 'regular':
        # There is something at the manifest's path and it is not a file this database
        # may execute -- a link out of the tree, or not a file at all. Refused, and
        # unanswered rather than unquarantined: the refusal is not a clean bill of health.
        _warn_unsafe_manifest(family_path, 'it is not a regular file')
        answer = (None, False)
    else:
        quarantine, parsed = _load_with_identity(label, family_path)
        if quarantine is None:
            # The manifest was there when this function looked and gone when the loader
            # did. Two different databases produce that, neither of them benign: a
            # concurrent edit, or a path that stopped resolving between the two reads.
            # Recording it as "carries no manifest" would turn a race into a clean bill
            # of health -- and cache it. It is an unanswered question instead, and is
            # deliberately NOT cached, because the next call may find the file back.
            logging.warning(
                'Quarantine manifest for family %s disappeared while it was being read '
                '(%s). Treating the family as unanswered rather than unquarantined; if '
                'the database is being edited while RMG runs, restart the run.',
                label, os.path.join(family_path, QUARANTINE_FILENAME))
            return None, False
        answer = (quarantine, True)
        # Cache under the identity of the file that was PARSED, not the one that was
        # stat'd above. The two are separate resolutions of one name, so a rename between
        # them would otherwise store B's quarantine under A's signature -- and restoring A
        # would then be a cache HIT returning B's answer, bypassing A's criterion without
        # reading A at all. Re-keying here is what makes the check and the thing used the
        # same object; it costs nothing, because the fstat was already taken to confirm
        # the file was regular.
        signature = (directory_exists, parsed, kind)
    _DISK_QUARANTINE_CACHE[key] = (signature, answer)
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
    key = (tuple(labels), type(kinetics).__name__)
    if key in _UNANSWERED_WARNED:
        # Asked before the suppression test on purpose: that test may now scan the
        # families directory rather than trust a cached negative, and this is what keeps
        # the scan to once per distinct question instead of once per edge reaction per
        # iteration.
        return
    if not _database_has_any_quarantine():
        # No quarantine anywhere in this database: nothing is being skipped, and warning
        # here would put a line into every ordinary run that loads a foreign seed.
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
