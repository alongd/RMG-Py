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

The governing ruling has two halves, and this file tests both.

**The gate must hard-fail, naming five things.** A run that admits a quarantined rate
must stop with an error naming the family, the rule or training-entry provenance, the
generated reaction, the kinetics class, and the reason. Not a warning, not a filtered
reaction, not a zeroed rate.

**Seven ways of making the run *appear* to succeed are forbidden**, each because it is
silent: removing the reaction after generation, substituting an average rule,
evaluating at ``potential = 0``, marking the reaction irreversible and continuing,
deriving a reverse rate, seeding a product or cation to bypass the channel, and
replacing it with a generic collision rate. There is one test per bullet below,
grouped in :class:`TestTheSevenForbiddenSilentBehaviours`.

The two halves pull in opposite directions and the third group,
:class:`TestTheDoorLeftOpen`, is where that tension is pinned: the ruling *also*
requires the family and its data to be preserved and to stay relevant to an
electrode/electrolyte model, so a gate that made these rules unusable everywhere would
overshoot. Database load, reaction generation, and the rate law itself are all
deliberately left working.

Most tests here build a synthetic quarantined family, so they run against any database.
The ones marked ``database`` check the real ``Cation_R_Recombination`` manifest and skip
with a loud reason on a database that predates it.

**Two kinds of red state, and the difference matters.** A campaign round is accepted on a
test shown failing before the repair and passing after. A test that fails at the base with
an `ImportError` or an `AttributeError` for a name the repair *adds* is not evidence about
behaviour -- it pins the shape of the repair, which is worth having and is not the same
thing. Those are labelled **structural** in each test's docstring and counted separately
in the round's findings; every finding also carries at least one **behavioural** test,
which imports only names the base already had, so its red state is the defect itself.
Round 107 relabelled round 105's four structural reds in place, at the manager's request,
rather than leaving the distinction in a findings document only.
"""

import errno
import inspect
import logging
import os
import pickle
import re
import subprocess
import sys
import textwrap

import pytest

from rmgpy import settings
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import KineticsLibrary, LibraryReaction
from rmgpy.data.kinetics.quarantine import (
    QUARANTINE_FILENAME,
    KineticsQuarantine,
    authoring_families,
    authoring_family,
    check_quarantine,
    describe_provenance,
    get_quarantine,
    load_family_quarantine,
    resolve_quarantine,
)
from rmgpy.exceptions import DatabaseError, QuarantinedKineticsError
from rmgpy.kinetics.arrhenius import Arrhenius, Marcus
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species, TransitionState

#: The family whose data is quarantined in the shipped database.
REAL_FAMILY = "Cation_R_Recombination"

#: What the shipped manifest declares, so the tests below fail loudly if the string
#: the ruling specified is ever quietly reworded.
REAL_STATE = "QUARANTINED FOR QUANTITATIVE PLASMA USE"
REAL_REASON = "electrochemical reference/domain unavailable"

#: The affected set as counted by the ruling. This is a cross-check on the criterion,
#: not the source of truth -- every test that uses it enumerates from the database
#: first and compares, so a database change moves the count and says so.
EXPECTED_RULES = 7
EXPECTED_TRAINING = 5

MANIFEST = """
name = "Fake_Quarantined_Family/quarantine"
state = "QUARANTINED FOR TESTING"
appliesToKineticsClass = "Marcus"
reason = "a reason that must reach the error message"
shortDesc = "short"
longDesc = "long"
"""


def make_marcus():
    """A Marcus model with parameters in the range the real family uses."""
    return Marcus(
        A=(1.73e06, "m^3/(mol*s)"),
        n=2,
        lmbd_i_coefs=[21824.5, -0.0341626, -0.0013254, 4.92966e-07],
        beta=(1.2e10, "1/m"),
        wr=(0, "kJ/mol"),
        wp=(0, "kJ/mol"),
        lmbd_o=(0, "J/mol"),
        comment="Estimated from node Root_2R->C",
    )


def make_reaction(family="Fake_Quarantined_Family"):
    """A two-to-one template reaction shaped like the one the real family generates."""
    reaction = TemplateReaction(
        reactants=[Species(label="Lip", molecule=[Molecule(smiles="[Li+]")]),
                   Species(label="CH3", molecule=[Molecule(smiles="[CH3]")])],
        products=[Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")])],
        family=family,
        reversible=True,
    )
    reaction.template = ["Root_2R->C"]
    return reaction


def write_manifest(directory, body=MANIFEST):
    """Write a manifest into `directory` and return the family path."""
    with open(os.path.join(directory, QUARANTINE_FILENAME), "w") as f:
        f.write(body)
    return directory


@pytest.fixture
def quarantine(tmp_path):
    """A loaded synthetic quarantine, from a real manifest file on disk."""
    return load_family_quarantine("Fake_Quarantined_Family", write_manifest(str(tmp_path)))


@pytest.fixture
def registered(monkeypatch, quarantine, tmp_path_factory):
    """
    Register the synthetic quarantined family in the kinetics database singleton.

    `check_quarantine` resolves a reaction's family label the same way the rest of RMG
    does, so this is what makes the gate see the synthetic family.

    The database directory is replaced too, with an empty one of its own. Enumeration
    reads the loaded families AND ``settings['database.directory']/kinetics/families`` on
    disk; replacing only the first left the second pointing at whatever database the
    checkout's rmgrc names, so a manifest merged into that database (an Arrhenius one, on
    2026-09-23) turned `test_an_ordinary_library_rate_says_nothing` red with no change to
    RMG-Py (round 113). A test that means to read a real database sets it itself.
    """
    database_directory = tmp_path_factory.mktemp("isolated_database")
    (database_directory / "kinetics" / "families").mkdir(parents=True)
    monkeypatch.setitem(settings, "database.directory", str(database_directory))

    class _Family:
        label = "Fake_Quarantined_Family"

        def __init__(self, q):
            self.quarantine = q

    class _Kinetics:
        def __init__(self, families):
            self.families = families

    class _Database:
        def __init__(self, families):
            self.kinetics = _Kinetics(families)

    import rmgpy.data.rmg

    families = {
        "Fake_Quarantined_Family": _Family(quarantine),
        "Ordinary_Family": _Family(None),
    }
    monkeypatch.setattr(rmgpy.data.rmg, "database", _Database(families), raising=False)
    return quarantine


class _StubModel(CoreEdgeReactionModel):
    """
    A reaction model whose kinetics estimation is stubbed out.

    `apply_kinetics_to_reaction` is the gate's primary site and the only caller of
    `generate_kinetics`; stubbing the latter isolates the gate from the database
    without touching the code path under test.
    """

    def __init__(self, result):
        super().__init__()
        self._result = result

    def generate_kinetics(self, reaction):
        return self._result


class TestTheManifest:
    """The marker itself: what it declares, and what it refuses to declare."""

    def test_a_manifest_loads_its_fields(self, quarantine):
        assert isinstance(quarantine, KineticsQuarantine)
        assert quarantine.state == "QUARANTINED FOR TESTING"
        assert quarantine.reason == "a reason that must reach the error message"
        assert quarantine.kinetics_class is Marcus
        assert quarantine.path.endswith(QUARANTINE_FILENAME)

    def test_a_family_without_a_manifest_is_not_quarantined(self, tmp_path):
        """The ordinary case, and the reason the gate costs nothing for real chemistry."""
        assert load_family_quarantine("Ordinary_Family", str(tmp_path)) is None

    def test_an_unknown_kinetics_class_raises_rather_than_quarantining_nothing(self, tmp_path):
        """
        A manifest that resolves to no class would report a quarantine while gating
        nothing -- strictly worse than having no manifest, because it reads as
        protection.
        """
        path = write_manifest(str(tmp_path), MANIFEST.replace('"Marcus"', '"Marcuss"'))
        with pytest.raises(DatabaseError) as exc:
            load_family_quarantine("Fake_Quarantined_Family", path)
        assert "Marcuss" in str(exc.value)

    def test_a_declared_engine_requirement_is_honoured(self, tmp_path):
        """
        A manifest may declare the engine capability its refusal depends on. If nothing
        reads that declaration it is a comment, not a pin: a database could state a
        requirement, be loaded by an engine that does not meet it, and gate nothing while
        the file still reads as protection.
        """
        body = MANIFEST + 'requiresEngineModule = "rmgpy.data.kinetics.no_such_module"\n'
        path = write_manifest(str(tmp_path), body)
        with pytest.raises(DatabaseError) as exc:
            load_family_quarantine("Fake_Quarantined_Family", path)
        assert "no_such_module" in str(exc.value)
        assert "Fake_Quarantined_Family" in str(exc.value)

    def test_a_declared_engine_symbol_is_honoured(self, tmp_path):
        """
        The drift that actually happens: the module survives a refactor and the function
        the gate depends on is renamed out from under the manifest.
        """
        body = (MANIFEST
                + 'requiresEngineModule = "rmgpy.data.kinetics.quarantine"\n'
                + 'requiresEngineSymbol = "a_function_this_engine_does_not_have"\n')
        path = write_manifest(str(tmp_path), body)
        with pytest.raises(DatabaseError) as exc:
            load_family_quarantine("Fake_Quarantined_Family", path)
        assert "a_function_this_engine_does_not_have" in str(exc.value)

    def test_a_satisfied_engine_requirement_loads(self, tmp_path):
        """The positive control: the check must not refuse a requirement that IS met."""
        body = (MANIFEST
                + 'requiresEngineModule = "rmgpy.data.kinetics.quarantine"\n'
                + 'requiresEngineSymbol = "check_quarantine"\n')
        path = write_manifest(str(tmp_path), body)
        assert load_family_quarantine("Fake_Quarantined_Family", path) is not None

    def test_the_commit_field_is_not_declarable(self, tmp_path):
        """
        An installed engine has no reliable commit to compare against, so a commit check
        would pass on every checkout: a check that cannot fail.

        Round 80 answered that by documenting the field as provenance and leaving it
        readable. That was the same defect one layer up -- a field still *named*
        ``requiresEngineCommit`` reads as a guarantee to anyone who greps the name rather
        than the comment beside it. The field is now refused outright, and the error says
        where to put the commit instead.
        """
        body = (MANIFEST
                + 'requiresEngineModule = "rmgpy.data.kinetics.quarantine"\n'
                + 'requiresEngineSymbol = "check_quarantine"\n'
                + 'requiresEngineCommit = "not-a-commit-that-exists-anywhere"\n')
        path = write_manifest(str(tmp_path), body)
        with pytest.raises(DatabaseError) as exc:
            load_family_quarantine("Fake_Quarantined_Family", path)
        assert "recordedEngineCommit" in str(exc.value)

    def test_the_recorded_commit_is_kept_as_provenance_and_reaches_a_reader(
            self, tmp_path, caplog):
        """
        The commit is still worth recording -- under a name that does not claim to pin,
        and read by something, because a field nothing consults is where this whole
        thread started.
        """
        body = (MANIFEST
                + 'requiresEngineModule = "rmgpy.data.kinetics.quarantine"\n'
                + 'requiresEngineSymbol = "check_quarantine"\n'
                + 'recordedEngineCommit = "541e6498f"\n')
        path = write_manifest(str(tmp_path), body)
        with caplog.at_level(logging.INFO):
            assert load_family_quarantine("Fake_Quarantined_Family", path) is not None
        assert "541e6498f" in caplog.text

    def test_a_symbol_that_is_not_callable_is_refused(self, tmp_path):
        """
        ``getattr(module, name) is not None`` accepts any attribute at all: a manifest
        naming ``math.pi`` as the gate that refuses a family used to load clean. A gate
        that cannot be called is not a gate.
        """
        body = (MANIFEST
                + 'requiresEngineModule = "math"\n'
                + 'requiresEngineSymbol = "pi"\n')
        path = write_manifest(str(tmp_path), body)
        with pytest.raises(DatabaseError) as exc:
            load_family_quarantine("Fake_Quarantined_Family", path)
        assert "not callable" in str(exc.value)

    def test_a_symbol_without_a_module_is_refused_not_ignored(self, tmp_path):
        """
        A symbol with nowhere to look it up asks for a check that cannot be performed.
        Skipping it silently leaves the manifest reading as pinned while pinning nothing,
        which is the failure mode of the whole field group.
        """
        body = MANIFEST + 'requiresEngineSymbol = "check_quarantine"\n'
        path = write_manifest(str(tmp_path), body)
        with pytest.raises(DatabaseError) as exc:
            load_family_quarantine("Fake_Quarantined_Family", path)
        assert "requiresEngineModule" in str(exc.value)

    def test_a_gate_that_exists_but_is_not_wired_in_is_refused(self, tmp_path):
        """
        Existence of the symbol shows the capability was written, not that anything calls
        it. ``os`` is a module that certainly exists and certainly does not bind the gate,
        so a manifest requiring the gate to be reached from there must be refused.
        """
        body = (MANIFEST
                + 'requiresEngineModule = "rmgpy.data.kinetics.quarantine"\n'
                + 'requiresEngineSymbol = "check_quarantine"\n'
                + 'requiresEngineCallSites = ("os",)\n')
        path = write_manifest(str(tmp_path), body)
        with pytest.raises(DatabaseError) as exc:
            load_family_quarantine("Fake_Quarantined_Family", path)
        assert "wired into" in str(exc.value)

    def test_a_call_site_that_imports_the_gate_and_never_calls_it_is_refused(
            self, tmp_path, monkeypatch):
        """
        Binding was the first version of this check and it is not the property the field
        is named for: delete every call in ``model.py`` and keep the import, and the
        module still binds the gate while running nothing. Round 87's ``math.pi`` finding
        one level up -- stricter about what the name points at, still silent about whether
        anything invokes it.
        """
        module_dir = tmp_path / "fake_site"
        module_dir.mkdir()
        (module_dir / "imports_but_never_calls.py").write_text(
            "from rmgpy.data.kinetics.quarantine import check_quarantine\n"
            "# every call site deleted; the import is all that is left\n")
        monkeypatch.syspath_prepend(str(module_dir))

        manifest_dir = tmp_path / "manifest"
        manifest_dir.mkdir()
        write_manifest(str(manifest_dir),
                       MANIFEST
                       + 'requiresEngineModule = "rmgpy.data.kinetics.quarantine"\n'
                       + 'requiresEngineSymbol = "check_quarantine"\n'
                       + 'requiresEngineCallSites = ("imports_but_never_calls",)\n')
        with pytest.raises(DatabaseError) as exc:
            load_family_quarantine("Fake_Quarantined_Family", str(manifest_dir))
        assert "never calls it" in str(exc.value)

    def test_the_real_call_site_satisfies_the_wiring_check(self, tmp_path):
        """
        The positive control for the test above, and the arrangement the shipped manifest
        declares: the reaction model imports the gate by name, so the binding is the same
        object.
        """
        body = (MANIFEST
                + 'requiresEngineModule = "rmgpy.data.kinetics.quarantine"\n'
                + 'requiresEngineSymbol = "check_quarantine"\n'
                + 'requiresEngineCallSites = ("rmgpy.rmg.model",)\n')
        path = write_manifest(str(tmp_path), body)
        assert load_family_quarantine("Fake_Quarantined_Family", path) is not None

    @pytest.mark.parametrize("field", ["state", "appliesToKineticsClass", "reason"])
    def test_a_manifest_missing_a_required_field_raises(self, tmp_path, field):
        body = "\n".join(line for line in MANIFEST.splitlines()
                         if not line.startswith(field + " "))
        path = write_manifest(str(tmp_path), body)
        with pytest.raises(DatabaseError) as exc:
            load_family_quarantine("Fake_Quarantined_Family", path)
        assert field in str(exc.value)

    def test_the_criterion_is_selective(self, quarantine):
        """
        Quarantining a family must not quarantine everything in it. If this ever
        passed vacuously -- by matching any kinetics at all -- every test below would
        pass for the wrong reason.
        """
        assert quarantine.applies_to(make_marcus())
        assert not quarantine.applies_to(Arrhenius(A=(1e13, "cm^3/(mol*s)"), n=0, Ea=(0, "kJ/mol")))
        assert not quarantine.applies_to(None)

    def test_the_criterion_is_evaluated_against_the_data_not_a_stored_list(self, quarantine):
        """
        The drift property. `affected_entries` reads the family it is handed, so an
        entry added after the manifest was written is covered without editing the
        manifest, and an entry refitted to another kinetics class leaves quarantine
        the same way. Nothing anywhere records *which* entries are affected.
        """

        class _Entry:
            def __init__(self, index, data):
                self.index, self.label, self.data = index, str(index), data

        class _Rules:
            entries = {"Root": [_Entry(1, make_marcus())],
                       "Root_2R->C": [_Entry(2, make_marcus())]}

        class _Depository:
            label = "Fake_Quarantined_Family/training"
            entries = {0: _Entry(0, make_marcus()),
                       1: _Entry(1, Arrhenius(A=(1e13, "cm^3/(mol*s)"), n=0, Ea=(0, "kJ/mol")))}

        class _Family:
            rules = _Rules()
            depositories = [_Depository()]

        family = _Family()
        affected = quarantine.affected_entries(family)
        assert len(affected["rules"]) == 2
        assert len(affected["training"]) == 1  # the Arrhenius entry is not quarantined

        # Add an entry the manifest has never heard of; it is covered immediately.
        _Rules.entries["Root_N-2R->C"] = [_Entry(3, make_marcus())]
        assert len(quarantine.affected_entries(family)["rules"]) == 3

        # Refit one to Arrhenius; it leaves quarantine immediately.
        _Rules.entries["Root"][0].data = Arrhenius(A=(1e13, "cm^3/(mol*s)"), n=0, Ea=(0, "kJ/mol"))
        assert len(quarantine.affected_entries(family)["rules"]) == 2

    def test_no_hard_coded_entry_list_exists_in_the_code(self):
        """
        The failure mode the ruling names by name. If a list of affected entries ever
        appears in the loader or the gate, it can drift out of step with the database,
        and this test is what says so.
        """
        import rmgpy.data.kinetics.quarantine as module

        source = inspect.getsource(module)
        for label in ("Root_2R->C", "Root_N-2R->C", "NH2 + Li", "C2H5 + Li"):
            assert label not in source, f"{label!r} is hard-coded in {module.__file__}"
        assert REAL_FAMILY not in source, f"{REAL_FAMILY!r} is hard-coded in {module.__file__}"


class TestTheRefusalNamesFiveThings:
    """The ruling's five required fields, one assertion each."""

    def _refusal(self, registered, **kwargs):
        reaction = make_reaction()
        with pytest.raises(QuarantinedKineticsError) as exc:
            check_quarantine(reaction, stage="a test", kinetics=make_marcus(), **kwargs)
        return str(exc.value)

    def test_it_names_the_family(self, registered):
        assert "Fake_Quarantined_Family" in self._refusal(registered)

    def test_it_names_the_provenance(self, registered):
        assert "Root_2R->C" in self._refusal(registered)

    def test_it_names_the_generated_reaction(self, registered):
        message = self._refusal(registered)
        assert "Lip" in message and "CH3Li" in message

    def test_it_names_the_kinetics_class(self, registered):
        assert "Marcus" in self._refusal(registered)

    def test_it_names_the_reason(self, registered):
        assert "a reason that must reach the error message" in self._refusal(registered)

    def test_it_names_the_campaign_state_and_the_manifest(self, registered):
        message = self._refusal(registered)
        assert "QUARANTINED FOR TESTING" in message
        assert QUARANTINE_FILENAME in message

    def test_provenance_is_never_blank(self, registered):
        """
        The provenance field has three sources of decreasing quality and no fourth.
        An averaged rate-rule estimate supplies no Entry at all, so if the fallbacks
        did not hold, this is the case where the refusal would name four fields.
        """
        reaction = make_reaction()
        reaction.template = None
        kinetics = make_marcus()
        kinetics.comment = ""
        assert describe_provenance(reaction, kinetics=kinetics) == "unrecorded"
        assert describe_provenance(reaction, kinetics=make_marcus())  # comment fallback
        assert describe_provenance(make_reaction(), kinetics=kinetics)  # template fallback


class TestTheGateFires:
    """Where the refusal happens, and where it deliberately does not."""

    def test_at_kinetics_estimation_the_primary_site(self, registered):
        model = _StubModel((make_marcus(), "rate rules", None, True))
        reaction = make_reaction()
        with pytest.raises(QuarantinedKineticsError):
            model.apply_kinetics_to_reaction(reaction)

    def test_at_core_admission_the_backstop(self, registered):
        model = CoreEdgeReactionModel()
        reaction = make_reaction()
        reaction.kinetics = make_marcus()
        with pytest.raises(QuarantinedKineticsError):
            model.add_reaction_to_core(reaction)
        assert reaction not in model.core.reactions

    def test_at_edge_admission_the_backstop(self, registered):
        """
        The edge is not a holding pen: edge fluxes decide what is promoted, so a rate
        of 1e-226 sitting there has already declared a real channel unimportant. The
        preflight deck this ticket came from put the Marcus reaction in the edge and
        never in the core, so a gate at core admission alone would not have fired.
        """
        model = CoreEdgeReactionModel()
        reaction = make_reaction()
        reaction.kinetics = make_marcus()
        with pytest.raises(QuarantinedKineticsError):
            model.add_reaction_to_edge(reaction)
        assert reaction not in model.edge.reactions

    def test_ordinary_chemistry_is_untouched(self, registered):
        """
        The negative control, in unit form: a family with no manifest is admitted
        normally even when its kinetics happen to be the quarantined class.
        """
        model = CoreEdgeReactionModel()
        reaction = make_reaction(family="Ordinary_Family")
        reaction.kinetics = make_marcus()
        model.add_reaction_to_core(reaction)
        model.add_reaction_to_edge(reaction)
        assert reaction in model.core.reactions

    def test_a_quarantined_family_still_admits_its_unquarantined_kinetics(self, registered):
        """The quarantine is a criterion, not a blanket ban on the family."""
        model = CoreEdgeReactionModel()
        reaction = make_reaction()
        reaction.kinetics = Arrhenius(A=(1e13, "cm^3/(mol*s)"), n=0, Ea=(0, "kJ/mol"))
        model.add_reaction_to_core(reaction)
        assert reaction in model.core.reactions

    def test_an_unloaded_database_does_not_crash_the_gate(self, monkeypatch):
        """
        `get_quarantine` is called on every reaction admission, including in unit
        tests and API use where no database has been loaded. It must return None, not
        raise -- a gate that breaks unrelated code would get switched off.
        """
        import rmgpy.data.rmg

        monkeypatch.setattr(rmgpy.data.rmg, "database", None, raising=False)
        assert get_quarantine("anything") is None
        check_quarantine(make_reaction(), stage="a test", kinetics=make_marcus())


def _clear_gate_caches():
    """
    Reset the gate's memoised state between tests.

    Each of these exists so a message is emitted once rather than once per edge reaction
    per iteration, or so a disk lookup is not repeated. Left uncleared, one test's warning
    silences another's and the silence reads as a pass.

    Deliberately tolerant of a cache that does not exist: this file is run against older
    engines to show new tests failing, and a fixture that raises `AttributeError` during
    setup turns an assertion failure into a setup error, which proves much less.
    """
    from rmgpy.data.kinetics import quarantine as module

    for name in ("_UNATTRIBUTED_WARNED", "_UNENUMERABLE_WARNED", "_UNANSWERED_WARNED",
                 "_DISK_QUARANTINE_CACHE",
                 "_DISK_ANY_QUARANTINE_CACHE", "_UNSAFE_LABELS_WARNED",
                 "_UNSAFE_MANIFESTS_WARNED", "_LEGACY_CALL_SITES_WARNED"):
        cache = getattr(module, name, None)
        if cache is not None:
            cache.clear()


def _register_families(monkeypatch, quarantines):
    """
    Register `quarantines` (label -> KineticsQuarantine or None) as the loaded families.

    Deliberately separate from the `registered` fixture so a test can register families
    that do NOT include the one a reaction names -- which is the case
    `add_seed_mechanism_to_core` creates routinely.
    """

    class _Family(object):
        def __init__(self, label, q):
            self.label = label
            self.quarantine = q

    class _Kinetics(object):
        def __init__(self, families):
            self.families = families

    class _Database(object):
        def __init__(self, families):
            self.kinetics = _Kinetics(families)

    import rmgpy.data.rmg

    families = {label: _Family(label, q) for label, q in quarantines.items()}
    monkeypatch.setattr(rmgpy.data.rmg, "database", _Database(families), raising=False)
    return families


def _deny_access(monkeypatch, directory):
    """
    Make `directory` answer the enumeration's system calls as a mode-000 directory does:
    it can itself be stat'd, but it cannot be listed and nothing below it can be reached.

    Injected by patch, not by ``chmod``. ``chmod 000`` does not stop root, so a permission
    test built on it has to skip under root -- and a test that can skip itself green
    proves nothing where it skips (round 113). Returns a callable that lifts the denial.
    """
    directory = os.path.abspath(str(directory))
    state = {"denied": True}

    def blocked(path, inside_only):
        if not state["denied"] or isinstance(path, int):
            return False
        path = os.path.abspath(os.fsdecode(path))
        return path.startswith(directory + os.sep) or (not inside_only and path == directory)

    def deny(name, inside_only):
        real = getattr(os, name)

        def denied(path=".", *args, **kwargs):
            if kwargs.get("dir_fd") is None and blocked(path, inside_only):
                raise PermissionError(errno.EACCES, os.strerror(errno.EACCES), path)
            return real(path, *args, **kwargs)

        monkeypatch.setattr(os, name, denied)

    # Not `os.open`: `_read_manifest` requires `os.open in os.supports_dir_fd`, which a
    # wrapper fails, and every path reaches an `lstat` before it opens anything.
    for name in ("stat", "lstat"):
        deny(name, True)
    for name in ("listdir", "scandir"):
        deny(name, False)
    return lambda: state.update(denied=False)


def make_library_reaction(library="copied_seed", comment="", long_desc=None):
    """A library reaction shaped like one loaded from a seed mechanism."""
    reaction = LibraryReaction(
        reactants=[Species(label="Lip", molecule=[Molecule(smiles="[Li+]")]),
                   Species(label="CH3", molecule=[Molecule(smiles="[CH3]")])],
        products=[Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")])],
        library=library,
        kinetics=make_marcus() if comment is None else _marcus_with(comment),
        reversible=True,
    )
    if long_desc is not None:
        reaction.entry = Entry(long_desc=long_desc)
    return reaction


def _marcus_with(comment):
    kinetics = make_marcus()
    kinetics.comment = comment
    return kinetics


#: The provenance RMG writes for an estimated rate, in the shape it survives into a
#: library entry's longDesc.
ESTIMATED_COMMENT = ("Estimated using template [Root_2R->C] for rate rule [Root_2R->C]\n"
                     "Euclidian distance = 0\n"
                     "family: Fake_Quarantined_Family")


class TestProvenanceNotTheFamilySlot:
    """
    A quarantine is about the RATE, scoped by the family that authored it -- so the gate
    must resolve *authorship*, and ``reaction.family`` is not authorship.

    ``LibraryReaction.__init__`` assigns ``self.family = library``: the same attribute
    holds a family label on one wrapper and a library label on another. Keying a gate on
    a slot a wrapper repurposes fails in both directions at once, and both are tested
    here -- the quarantined rate that gets in, and the innocent library that gets refused.
    """

    @pytest.fixture(autouse=True)
    def _fresh_warning_cache(self):
        """
        The unattributable-rate warning is emitted once per (library, kinetics class) so
        it does not repeat per edge reaction per iteration. Clear it, or one test's
        warning silences another's and the silence reads as a pass.
        """
        _clear_gate_caches()
        yield
        _clear_gate_caches()

    def test_a_quarantined_rate_copied_into_a_library_is_still_refused(self, registered):
        """
        The HIGH. Copy the rate into a seed mechanism and its ``.family`` becomes
        ``"copied_seed"``, which resolves to no quarantine at all.
        """
        reaction = make_library_reaction(comment=ESTIMATED_COMMENT)
        with pytest.raises(QuarantinedKineticsError) as exc:
            check_quarantine(reaction, stage="a test", kinetics=reaction.kinetics)
        assert "Fake_Quarantined_Family" in str(exc.value)

    def test_the_provenance_may_live_in_the_entry_instead_of_the_comment(self, registered):
        """
        The library writer puts the estimator's comment in the entry's longDesc.

        NOTE the shape this builds: it assigns ``reaction.entry`` by hand. That is the
        intended design, and for one round it was a design production did not implement --
        ``get_library_reactions`` constructed every ``LibraryReaction`` without its entry,
        so this test passed against an object shape that never occurred. It is kept as a
        unit test of the lookup, and
        :meth:`TestTheProvenanceArrivesFromTheRealLoader.test_the_entry_travels_with_the_reaction`
        is what makes it mean anything.
        """
        reaction = make_library_reaction(comment="", long_desc=ESTIMATED_COMMENT)
        with pytest.raises(QuarantinedKineticsError):
            check_quarantine(reaction, stage="a test", kinetics=reaction.kinetics)

    def test_every_declared_family_is_consulted_not_just_the_first(self, registered):
        """
        Provenance is free text. Taking the first ``family:`` line and stopping means one
        prepended line shadows the genuine one, and a bypass costing an attacker a single
        comment line is not a bound worth having. All declared labels are consulted, so an
        added line can only widen what is checked.
        """
        shadowed = make_library_reaction(
            comment="family: An_Innocent_Family\n" + ESTIMATED_COMMENT)
        assert authoring_families(shadowed) == ["An_Innocent_Family",
                                                "Fake_Quarantined_Family"]
        with pytest.raises(QuarantinedKineticsError):
            check_quarantine(shadowed, stage="a test", kinetics=shadowed.kinetics)

    def test_a_forged_label_is_bounded_by_the_criterion(self, registered):
        """
        The other edge of trusting free text, stated so the bound is visible: a forged or
        stale line can cause a false *refusal*, never a false admission, and only for a
        rate that already matches the manifest's kinetics criterion. A forged label on an
        Arrhenius rate does nothing.
        """
        forged = make_library_reaction(library="lib")
        forged.kinetics = Arrhenius(A=(1e13, "cm^3/(mol*s)"), n=0, Ea=(0, "kJ/mol"),
                                    comment=ESTIMATED_COMMENT)
        assert check_quarantine(forged, stage="a test",
                                kinetics=forged.kinetics) is None

    def test_a_library_named_like_a_quarantined_family_is_not_refused(self, registered):
        """
        The other direction, and the one a false-positive-hunting reviewer would find
        first: an unrelated library that happens to share the family's name carries no
        quarantined data, and refusing it would obstruct a legitimate mechanism.
        """
        reaction = make_library_reaction(library="Fake_Quarantined_Family", comment="")
        assert check_quarantine(reaction, stage="a test",
                                kinetics=reaction.kinetics) is None

    def test_authoring_family_reads_authorship_not_the_slot(self, registered):
        template = make_reaction()
        assert authoring_family(template) == "Fake_Quarantined_Family"

        copied = make_library_reaction(comment=ESTIMATED_COMMENT)
        assert copied.family == "copied_seed"
        assert authoring_family(copied) == "Fake_Quarantined_Family"

        anonymous = make_library_reaction(library="Fake_Quarantined_Family", comment="")
        assert anonymous.family == "Fake_Quarantined_Family"
        assert authoring_family(anonymous) is None

    def test_an_unattributable_rate_is_reported_and_not_refused(self, registered, caplog):
        """
        Where provenance is genuinely lost -- a hand-written entry with no comment -- the
        gate cannot tell a copied quarantined rate from an independent one of the same
        class. It says so once rather than guessing in either direction, because refusing
        on the kinetics class alone would ban legitimate electrochemistry.
        """
        reaction = make_library_reaction(library="some_hand_written_library", comment="")
        with caplog.at_level(logging.WARNING):
            assert check_quarantine(reaction, stage="a test",
                                    kinetics=reaction.kinetics) is None
        assert "some_hand_written_library" in caplog.text
        assert "NOT a refusal" in caplog.text

    def test_an_ordinary_library_rate_says_nothing(self, registered, caplog):
        """
        The control that keeps the warning from being noise: a library rate that does not
        match any live quarantine criterion must produce no output at all.
        """
        assert not os.listdir(os.path.join(settings["database.directory"], "kinetics",
                                           "families")), (
            "the fixture's database is not its own: this control must not depend on what "
            "the checked-out database happens to quarantine")
        reaction = make_library_reaction(library="primaryH2O2", comment="")
        reaction.kinetics = Arrhenius(A=(1e13, "cm^3/(mol*s)"), n=0, Ea=(0, "kJ/mol"))
        with caplog.at_level(logging.WARNING):
            assert check_quarantine(reaction, stage="a test",
                                    kinetics=reaction.kinetics) is None
        assert caplog.text == ""


class TestTheProvenanceArrivesFromTheRealLoader:
    """
    The key is right; these are the paths that have to deliver it.

    Every test here goes through a **production** constructor. A gate keyed on a field
    that production never fills is a gate that reads a design document, and the unit tests
    above cannot see the difference -- they build the field themselves.
    """

    @pytest.fixture(autouse=True)
    def _fresh_caches(self):
        _clear_gate_caches()
        yield
        _clear_gate_caches()

    @staticmethod
    def _library(label="copied_seed", long_desc=ESTIMATED_COMMENT, comment=""):
        """A real KineticsLibrary with one real Entry, loaded the way RMG loads one."""
        from rmgpy.reaction import Reaction

        library = KineticsLibrary(label=label)
        library.entries = {
            1: Entry(
                index=1,
                label="Lip + CH3 <=> CH3Li",
                item=Reaction(
                    reactants=[Species(label="Lip", molecule=[Molecule(smiles="[Li+]")]),
                               Species(label="CH3", molecule=[Molecule(smiles="[CH3]")])],
                    products=[Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")])],
                    reversible=True),
                data=_marcus_with(comment),
                long_desc=long_desc,
            )
        }
        return library

    def test_the_entry_travels_with_the_reaction(self, registered):
        """
        The one that could not fail before. `get_library_reactions` built every
        `LibraryReaction` without `entry=entry`, so the longDesc half of the lookup could
        never fire in production no matter how correct it was.
        """
        reaction = self._library().get_library_reactions()[0]

        assert reaction.entry is not None, (
            "the loader dropped the entry, so the provenance it carries is unreachable "
            "from the reaction and every longDesc lookup is dead code in production")
        assert authoring_family(reaction) == "Fake_Quarantined_Family"
        with pytest.raises(QuarantinedKineticsError):
            check_quarantine(reaction, stage="a test", kinetics=reaction.kinetics)

    def test_the_entry_survives_the_pickle_the_parallel_path_uses(self, registered):
        """
        `LibraryReaction.__reduce__` has to carry the entry too, or provenance is lost
        again the moment a reaction crosses into a worker process.
        """
        import pickle

        reaction = self._library().get_library_reactions()[0]
        restored = pickle.loads(pickle.dumps(reaction))
        assert restored.entry is not None
        assert authoring_family(restored) == "Fake_Quarantined_Family"

    def test_an_ordinary_library_still_loads_with_no_authorship(self, registered):
        """The control: a library with no provenance must still produce reactions."""
        library = self._library(label="primaryH2O2", long_desc="a hand-written entry")
        reaction = library.get_library_reactions()[0]
        assert reaction.entry is not None
        assert authoring_family(reaction) is None


class TestAMissingFamilyIsAnUnansweredQuestion:
    """
    ``None`` from "no authorship recorded" and ``None`` from "authorship recovered, family
    not found" are the same value arriving from opposite situations, and a gate that
    cannot tell them apart admits on a lookup miss.

    ``add_seed_mechanism_to_core`` makes the second case routine: it *converts* a reaction
    whose family is unavailable into a library reaction rather than loading the family.
    """

    @pytest.fixture(autouse=True)
    def _fresh_caches(self):
        _clear_gate_caches()
        yield
        _clear_gate_caches()

    @staticmethod
    def _reaction():
        return make_library_reaction(library="a_seed_from_another_chemistry",
                                     comment=ESTIMATED_COMMENT)

    def test_an_unloaded_family_is_answered_from_disk(self, monkeypatch, tmp_path, quarantine):
        """
        A manifest is a sidecar file; reading it needs no family object. So "the family
        exists in this database and merely was not loaded" is an answerable question, and
        answering it is what keeps a routine seed conversion from disabling the gate.
        """
        family_dir = tmp_path / "kinetics" / "families" / "Fake_Quarantined_Family"
        family_dir.mkdir(parents=True)
        write_manifest(str(family_dir))

        _register_families(monkeypatch, {"Some_Other_Loaded_Family": quarantine})
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))

        reaction = self._reaction()
        assert authoring_family(reaction) == "Fake_Quarantined_Family"
        with pytest.raises(QuarantinedKineticsError):
            check_quarantine(reaction, stage="a test", kinetics=reaction.kinetics)

    def test_a_family_in_the_database_without_a_manifest_is_a_clean_answer(
            self, monkeypatch, tmp_path, quarantine):
        """The negative control for the test above: present, unquarantined, admitted."""
        (tmp_path / "kinetics" / "families" / "Fake_Quarantined_Family").mkdir(parents=True)

        _register_families(monkeypatch, {"Some_Other_Loaded_Family": quarantine})
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))

        reaction = self._reaction()
        assert check_quarantine(reaction, stage="a test",
                                kinetics=reaction.kinetics) is None

    def test_an_unanswerable_question_is_reported_not_swallowed(
            self, monkeypatch, tmp_path, quarantine, caplog):
        """
        When the family is nowhere this run can reach, the answer is genuinely
        unavailable. It is admitted -- refusing would stop ordinary runs whose seeds name
        foreign families, and would break the promise that a database with no manifest
        behaves as before -- but it must not pass in silence, and the message must not be
        the unattributable one, because the authorship is right there.
        """
        _register_families(monkeypatch, {"Some_Other_Loaded_Family": quarantine})
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))

        reaction = self._reaction()
        with caplog.at_level(logging.WARNING):
            assert check_quarantine(reaction, stage="a test",
                                    kinetics=reaction.kinetics) is None
        assert "Fake_Quarantined_Family" in caplog.text
        assert "cannot know" in caplog.text
        assert "records no authoring family" not in caplog.text

    def test_a_database_with_no_quarantine_at_all_stays_silent(
            self, monkeypatch, tmp_path, caplog):
        """
        The bound on that warning. A database carrying no manifest anywhere must behave
        exactly as it did before -- otherwise every ordinary run that loads a foreign seed
        gains a line of output about a mechanism it does not use.
        """
        _register_families(monkeypatch, {"Some_Other_Loaded_Family": None})
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))

        with caplog.at_level(logging.WARNING):
            assert check_quarantine(self._reaction(), stage="a test",
                                    kinetics=make_marcus()) is None
        assert caplog.text == ""


class TestTheFamilySlotIsLeftAloneForItsOtherReader:
    """
    ``Reaction.family`` has a second reader that wants the OPPOSITE thing.

    ``rmgpy.electron_placement`` keys ``FAMILY_ELECTRON_PLACEMENT`` on
    ``Reaction.family`` and **depends** on ``LibraryReaction`` overwriting it:
    ``PlasmaElectronImpactIonization`` is a kinetics LIBRARY label sitting in that table
    deliberately, and the argon ionisation channel this whole campaign runs on resolves
    through it. The quarantine gate must not see that label as a family; the placement
    resolver must.

    Both are served only because :func:`authoring_family` leaves the slot alone and reads
    provenance beside it. These tests exist so that a later "cleanup" -- normalising
    ``family``, stopping ``LibraryReaction`` overwriting it, or promoting a separate
    authored-family attribute to authoritative -- breaks here instead of breaking argon
    ionisation placement silently, which nothing in this ticket would otherwise catch.
    """

    @staticmethod
    def _ionisation_library():
        """A real library labelled like the one the argon channel resolves through."""
        from rmgpy.reaction import Reaction

        library = KineticsLibrary(label="PlasmaElectronImpactIonization")
        library.entries = {
            1: Entry(
                index=1,
                label="Ar <=> Arp",
                item=Reaction(
                    reactants=[Species(label="Ar", molecule=[Molecule(smiles="[Ar]")])],
                    products=[Species(label="Arp",
                                      molecule=[Molecule().from_adjacency_list(
                                          "1 Ar u1 p3 c+1")])],
                    electrons=1,
                    reversible=False,
                ),
                data=Arrhenius(A=(1.254444e3, "m^3/(mol*s)"), n=0.0, Ea=(0.0, "kJ/mol")),
                long_desc="",
            )
        }
        return library

    def test_the_library_label_still_lands_in_the_family_slot(self):
        """
        The contract electron placement depends on, asserted through the loader this
        round changed. Adding `entry=entry` must not disturb what `family` holds.
        """
        reaction = self._ionisation_library().get_library_reactions()[0]
        assert reaction.library == "PlasmaElectronImpactIonization"
        assert reaction.family == "PlasmaElectronImpactIonization", (
            "electron placement is keyed on Reaction.family and this library label is a "
            "deliberate key in FAMILY_ELECTRON_PLACEMENT; changing what the slot holds "
            "breaks argon ionisation placement")
        assert reaction.entry is not None, "and the entry must arrive as well"

    def test_placement_still_resolves_for_the_argon_ionisation_channel(self):
        """
        End to end through the real resolver: the library label must still find its
        ``(1, 2)`` declaration and produce a balanced view.
        """
        from rmgpy.electron_placement import (
            FAMILY_ELECTRON_PLACEMENT,
            resolve_electron_placement,
        )

        assert FAMILY_ELECTRON_PLACEMENT["PlasmaElectronImpactIonization"] == (1, 2)

        reaction = self._ionisation_library().get_library_reactions()[0]
        electron = Species(label="e-").from_adjacency_list("1 e u1 p0 c-1")
        view = resolve_electron_placement(reaction, [electron] + reaction.reactants
                                          + reaction.products)

        assert sum(1 for spc in view.reactants if spc.is_electron()) == 1
        assert sum(1 for spc in view.products if spc.is_electron()) == 2
        assert view.electrons == 0
        # the canonical reaction is never mutated
        assert reaction.electrons == 1
        assert not any(spc.is_electron() for spc in reaction.reactants)

    def test_the_two_readers_disagree_about_the_same_object_and_both_are_right(self):
        """
        The design, pinned on one object. Placement reads the slot and gets the library;
        the quarantine gate reads provenance and gets the authoring family. Neither
        answer is a normalisation of the other, and making them agree would break one.
        """
        library = self._ionisation_library()
        library.entries[1].long_desc = "family: Plasma_Electron_Impact_Ionization"
        reaction = library.get_library_reactions()[0]

        assert reaction.family == "PlasmaElectronImpactIonization"
        assert authoring_family(reaction) == "Plasma_Electron_Impact_Ionization"
        assert reaction.family != authoring_family(reaction)

    def test_nothing_in_the_gate_writes_to_the_family_slot(self):
        """
        A static guard on the whole quarantine module: it may read `family`, never assign
        it. An assignment here is how the slot would drift under the resolver.
        """
        import ast

        from rmgpy.data.kinetics import quarantine as module

        tree = ast.parse(inspect.getsource(module))
        writes = [node for node in ast.walk(tree)
                  if isinstance(node, (ast.Assign, ast.AugAssign))
                  for target in (node.targets if isinstance(node, ast.Assign)
                                 else [node.target])
                  if isinstance(target, ast.Attribute) and target.attr == "family"]
        assert not writes, (
            "the quarantine module assigns to a .family attribute; it must only ever "
            "read the slot, because electron placement depends on what LibraryReaction "
            "puts there")


class TestEveryAdmissionPathIsGated:
    """
    The quarantine module's docstring enumerates every path by which a rate reaches the
    model and names the gate covering each. These tests are what make that table a claim
    rather than a comment.
    """

    def test_the_pressure_dependent_network_path_is_gated(self, registered):
        """
        A pdep path reaction goes to a network *instead of* core or edge, and
        ``generate_kinetics=False`` skips the estimation gate, so neither backstop sees
        it. Before this it reached a network unchecked and left as a k(T,P) fit with its
        provenance averaged away.
        """
        model = CoreEdgeReactionModel()
        reaction = make_reaction()
        reaction.kinetics = make_marcus()
        with pytest.raises(QuarantinedKineticsError):
            model.add_reaction_to_unimolecular_networks(
                reaction, new_species=reaction.reactants[0])

    def test_the_enumeration_in_the_docstring_names_every_gate_call(self):
        """
        Each gate named in the module docstring's table must be a real call in
        ``rmgpy/rmg/model.py``. A table that drifts from the code is worse than none.
        """
        import ast

        from rmgpy.data.kinetics import quarantine as module

        tree = ast.parse(inspect.getsource(CoreEdgeReactionModel))
        gated = {node.name for node in ast.walk(tree)
                 if isinstance(node, ast.FunctionDef)
                 and any(isinstance(call, ast.Call)
                         and getattr(call.func, "id", None) == "check_quarantine"
                         for call in ast.walk(node))}

        assert gated == {"apply_kinetics_to_reaction", "add_reaction_to_core",
                         "add_reaction_to_edge", "add_reaction_to_unimolecular_networks"}, (
            "the set of gated methods moved; update the table in the quarantine module "
            "docstring in the same commit, or it stops being checkable")
        for name in gated:
            assert name in module.__doc__, (
                "%s calls the gate but the docstring's admission-path table does not "
                "name it" % name)


class TestTheSevenForbiddenSilentBehaviours:
    """
    One test per bullet of the ruling's "Do not silently" list.

    Every one of these is a way to make a run *appear* to succeed. They are forbidden
    because they are silent, not because the outcome is necessarily wrong, so what is
    asserted throughout is that the exception escapes and that the reaction is left
    exactly as it was found.
    """

    def test_1_the_reaction_is_not_removed_after_generation(self, registered):
        """
        Forbidden: remove the reaction after generation.

        The gate raises; it never returns a sentinel that a caller could read as
        "skip this one", and the exception is not swallowed on the way out of the
        model. The reaction object survives untouched -- it is evidence, not litter.
        """
        model = _StubModel((make_marcus(), "rate rules", None, True))
        reaction = make_reaction()
        before = (list(reaction.reactants), list(reaction.products))

        with pytest.raises(QuarantinedKineticsError):
            model.apply_kinetics_to_reaction(reaction)

        assert (list(reaction.reactants), list(reaction.products)) == before
        assert model.core.reactions == [] and model.edge.reactions == []
        assert model.new_reaction_list == []
        # And the gate itself has no non-raising exit that reports a refusal.
        assert check_quarantine(make_reaction(family="Ordinary_Family"), stage="t") is None

    def test_2_no_average_rule_is_substituted(self, registered):
        """
        Forbidden: substitute an average rule.

        RMG's own averaged rate-rule estimate is the tempting substitute, because it
        arrives through the same call and looks like an ordinary result: an averaged
        estimate returns ``entry=None``, so a gate keyed on having exact provenance
        would let precisely this case through. It refuses identically.
        """
        model = _StubModel((make_marcus(), "rate rules", None, True))
        with pytest.raises(QuarantinedKineticsError) as exc:
            model.apply_kinetics_to_reaction(make_reaction())
        assert "Marcus" in str(exc.value)

        exact_entry = type("E", (), {"index": 3, "label": "Root_2R->C", "rank": 11})()
        model = _StubModel((make_marcus(), "rate rules", exact_entry, True))
        with pytest.raises(QuarantinedKineticsError):
            model.apply_kinetics_to_reaction(make_reaction())

    def test_3_it_is_never_evaluated_at_potential_zero(self, registered):
        """
        Forbidden: evaluate it with ``potential = 0``.

        The refusal happens before the kinetics are bound to the reaction, so there
        is nothing for any consumer to evaluate at any potential: `kinetics` is still
        None afterwards. The gate itself takes no potential and reads none, which is
        why no value of one can change its answer.
        """
        model = _StubModel((make_marcus(), "rate rules", None, True))
        reaction = make_reaction()
        assert reaction.kinetics is None

        with pytest.raises(QuarantinedKineticsError):
            model.apply_kinetics_to_reaction(reaction)

        assert reaction.kinetics is None
        assert "potential" not in inspect.signature(check_quarantine).parameters

    def test_4_it_is_not_marked_irreversible(self, registered):
        """
        Forbidden: mark it irreversible and continue.

        Two claims: `reversible` is not mutated, and there is no "and continue" --
        the exception leaves the method.
        """
        model = _StubModel((make_marcus(), "rate rules", None, True))
        reaction = make_reaction()
        assert reaction.reversible is True

        with pytest.raises(QuarantinedKineticsError):
            model.apply_kinetics_to_reaction(reaction)

        assert reaction.reversible is True

    def test_5_no_reverse_rate_is_derived(self, registered):
        """
        Forbidden: derive a reverse rate.

        This is the live case, not a hypothetical: the deck this ticket came from
        reaches ``Cation_R_Recombination`` through its auto-derived reverse template,
        so the kinetics match in the reverse direction and RMG's next act would be to
        flip the reaction and adopt the rate. The gate is placed before that flip, so
        reactants and products come out in the order they went in.
        """
        model = _StubModel((make_marcus(), "rate rules", None, False))  # is_forward=False
        reaction = make_reaction()
        reactants, products = list(reaction.reactants), list(reaction.products)

        with pytest.raises(QuarantinedKineticsError):
            model.apply_kinetics_to_reaction(reaction)

        assert list(reaction.reactants) == reactants
        assert list(reaction.products) == products
        assert reaction.kinetics is None

    def test_6_seeding_a_product_or_cation_does_not_bypass_it(self, registered):
        """
        Forbidden: seed a product or cation to bypass it.

        Seeding is how a run gets a quarantined reaction into the model without going
        through kinetics estimation -- a seed mechanism or reaction library carries
        its own kinetics, so `apply_kinetics_to_reaction` is never called for it.
        That is what the core and edge backstops are for, and this is the test that
        holds them in place.
        """
        for admit in ("add_reaction_to_core", "add_reaction_to_edge"):
            model = CoreEdgeReactionModel()
            reaction = make_reaction()
            reaction.kinetics = make_marcus()  # pre-attached, as from a seed mechanism
            with pytest.raises(QuarantinedKineticsError):
                getattr(model, admit)(reaction)
        assert model.core.reactions == [] and model.edge.reactions == []

    def test_7_no_generic_collision_rate_replaces_it(self, registered):
        """
        Forbidden: replace it with a generic collision rate.

        The gate has one non-raising outcome and it is ``None``. It cannot return a
        substitute because it returns no kinetics at all, and the refused reaction
        keeps the kinetics it had -- none.
        """
        model = _StubModel((make_marcus(), "rate rules", None, True))
        reaction = make_reaction()
        with pytest.raises(QuarantinedKineticsError):
            model.apply_kinetics_to_reaction(reaction)
        assert reaction.kinetics is None

        signature = inspect.signature(check_quarantine)
        assert signature.return_annotation is inspect.Signature.empty
        assert check_quarantine(make_reaction(family="Ordinary_Family"), stage="t") is None


@pytest.mark.database
class TestTheRealQuarantinedFamily:
    """
    The same guarantees against the family that actually ships the quarantined data.

    Skips loudly rather than failing on a database that predates the manifest, since
    the code worktree and the data worktree merge separately.
    """

    @classmethod
    def setup_class(cls):
        families_path = os.path.join(settings["database.directory"], "kinetics", "families")
        manifest = os.path.join(families_path, REAL_FAMILY, QUARANTINE_FILENAME)
        if not os.path.exists(manifest):
            pytest.skip(f"database at {settings['database.directory']} carries no quarantine "
                        f"manifest for {REAL_FAMILY} ({manifest} does not exist)")
        database = KineticsDatabase()
        database.load_families(path=families_path, families=[REAL_FAMILY])
        cls.family = database.families[REAL_FAMILY]

    def test_the_family_is_quarantined_with_the_ruling_s_words(self):
        assert self.family.quarantine is not None
        assert self.family.quarantine.state == REAL_STATE
        assert self.family.quarantine.reason == REAL_REASON
        assert self.family.quarantine.kinetics_class is Marcus

    def test_the_affected_set_is_enumerated_from_the_database(self):
        """
        Counts come from the data, and the ruling's numbers are the cross-check. If
        the database gains or loses a Marcus entry, this fails and names the new
        count rather than quietly disagreeing with the ledger.
        """
        affected = self.family.quarantine.affected_entries(self.family)
        assert len(affected["rules"]) == EXPECTED_RULES, \
            f"{len(affected['rules'])} rate rules affected, ruling recorded {EXPECTED_RULES}"
        assert len(affected["training"]) == EXPECTED_TRAINING, \
            f"{len(affected['training'])} training entries affected, " \
            f"ruling recorded {EXPECTED_TRAINING}"
        every = affected["rules"] + affected["training"]
        assert every, (
            "nothing was enumerated, and the class assertion below would then pass "
            "vacuously -- `all()` over an empty list is the shape this campaign keeps "
            "shipping as evidence")
        assert all(isinstance(e.data, Marcus) for e in every)

    def test_nothing_was_deleted(self):
        """
        Preservation is an explicit requirement: the entries are provenance evidence
        for the recovery track. Quarantine is a label plus a gate, never a removal.
        """
        all_rules = [entry for entries in self.family.rules.entries.values() for entry in entries]
        training = [d for d in self.family.depositories if d.label.endswith("/training")]
        assert len(all_rules) == EXPECTED_RULES
        assert training and len(training[0].entries) == EXPECTED_TRAINING

    def test_the_rates_themselves_were_not_touched(self):
        """
        No refit, no reinterpretation. Spot-check the parameters the ruling forbids
        adjusting on the node the preflight deck actually matches.
        """
        entry = self.family.rules.entries["Root_2R->C"][0]
        assert entry.data.A.value_si == pytest.approx(1.73e06)
        assert entry.data.n.value_si == pytest.approx(2)
        assert entry.data.lmbd_o.value_si == pytest.approx(0.0)
        assert entry.data.beta.value_si == pytest.approx(1.2e10)
        assert entry.data.wr.value_si == pytest.approx(0.0)
        assert entry.data.wp.value_si == pytest.approx(0.0)
        # The WHOLE polynomial, not its leading term. `lmbd_i_coefs` has four components
        # and a refit moves the later ones hardest -- checking [0] alone is a collapsed
        # assertion on a vector quantity, which passes for three of the four ways this
        # could be rewritten.
        coefficients = list(entry.data.lmbd_i_coefs.value_si)
        assert len(coefficients) == 4, \
            f"the reorganisation barrier polynomial has {len(coefficients)} coefficients, " \
            f"and the ruling was recorded against four"
        assert coefficients == pytest.approx(
            [51487.7, -0.166019, -0.00176034, 4.42738e-07], rel=1e-6)


@pytest.mark.database
class TestTheDoorLeftOpen:
    """
    What the quarantine deliberately does *not* break.

    The ruling preserves the family as potentially relevant to an electrode/electrolyte
    model, so the gate sits inside RMG's mechanism builder and nowhere else. Everything
    below is a consumer that keeps working on the same unmodified data.
    """

    @classmethod
    def setup_class(cls):
        families_path = os.path.join(settings["database.directory"], "kinetics", "families")
        if not os.path.exists(os.path.join(families_path, REAL_FAMILY, QUARANTINE_FILENAME)):
            pytest.skip(f"database at {settings['database.directory']} carries no quarantine "
                        f"manifest for {REAL_FAMILY}")
        database = KineticsDatabase()
        database.load_families(path=families_path, families=[REAL_FAMILY])
        cls.family = database.families[REAL_FAMILY]

    def test_the_database_still_loads(self):
        """A gate at database load would make the data unusable in every context."""
        assert self.family.rules.entries
        assert self.family.groups is not None

    def test_the_family_is_still_registered_and_still_generates_reactions(self):
        reactions = self.family.generate_reactions(
            [Molecule(smiles="[Li+]"), Molecule(smiles="[CH3]")])
        assert reactions, "the quarantined family stopped generating reactions"

    def test_the_rate_law_still_evaluates(self):
        """
        `Marcus.get_rate_coefficient` and `Reaction.get_rate_coefficient` are
        untouched. A consumer that supplies the missing electrochemical reference
        needs them, and refusing there would foreclose the recovery track this
        quarantine exists to keep open.
        """
        kinetics = self.family.rules.entries["Root_2R->C"][0].data
        assert kinetics.get_rate_coefficient(1000.0, -1.0e4) > 0.0


# ---------------------------------------------------------------------------------------
# Round 92
# ---------------------------------------------------------------------------------------


def _auto_generated_library(label, family_label, auto=True, comment=""):
    """
    A real `KineticsLibrary` whose one entry produces the auto-generated TEMPLATE shape.

    `get_library_reactions` builds three reaction shapes and this is the third: an entry
    whose longDesc says "rate rule" gets a `TemplateReaction` with the authoring family
    parsed out of that longDesc. The kinetics comment is left without a `family:` line on
    purpose, so the entry is the only carrier -- which is the case the other two shapes
    had covered since round 89 and this one did not.
    """
    library = KineticsLibrary(label=label, name=label)
    library.auto_generated = auto
    library.entries = {
        1: Entry(
            index=1,
            label="Lip + CH3 <=> CH3Li",
            item=Reaction(
                reactants=[Species(label="Lip", molecule=[Molecule(smiles="[Li+]")],
                                   reactive=False),
                           Species(label="CH3", molecule=[Molecule(smiles="[CH3]")],
                                   reactive=False)],
                products=[Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")],
                                  reactive=False)],
                reversible=False),
            data=_marcus_with(comment),
            long_desc=("Matched reaction 3 Lip + CH3 <=> CH3Li in {0}/rate rule [Root]\n"
                       "Euclidian distance = 0\n"
                       "family: {0}".format(family_label)),
        )
    }
    return library


#: A family that is quarantined in the database on disk and NOT loaded -- which is the
#: state that makes `CoreEdgeReactionModel` convert a template reaction in the first
#: place, and the state in which the gate has to answer from disk.
QUARANTINED_ON_DISK = "A_Family_Quarantined_On_Disk"


def _quarantined_on_disk(root, label=QUARANTINED_ON_DISK):
    """Lay out `root` as a database directory whose `label` family carries a manifest."""
    family = root / "kinetics" / "families" / label
    family.mkdir(parents=True)
    write_manifest(family)
    return family


def _model_database(monkeypatch, libraries, families):
    """
    The pieces `add_seed_mechanism_to_core` and `add_reaction_library_to_edge` read.

    Duck-typed, and the duck-typing is the point of the two tests below being run
    against the REAL methods rather than a copy of them: everything the methods touch on
    the way to the gate is here, and nothing else is stubbed.
    """

    class _Kinetics(object):
        def __init__(self):
            self.libraries = libraries
            self.families = families
            self.library_order = []

        def load_libraries(self, path=None, libraries=None):
            raise AssertionError("every library in this test is already known")

    class _Forbidden(object):
        def is_molecule_forbidden(self, molecule):
            return False

    class _Database(object):
        def __init__(self):
            self.kinetics = _Kinetics()
            self.forbidden_structures = _Forbidden()

    import rmgpy.data.rmg

    monkeypatch.setattr(rmgpy.data.rmg, "database", _Database(), raising=False)


class TestTheThirdShapeSurvivesTheConversion:
    """
    Round 92. `get_library_reactions` builds three reaction shapes; round 89 attached the
    authoring entry to the two `LibraryReaction` ones. The third is a `TemplateReaction`,
    and the two sites in `rmgpy/rmg/model.py` that rebuild one as a `LibraryReaction` when
    its family is not loaded carried neither the entry nor the parsed family across -- one
    statement after logging that family to the user. A rate the gate refused before the
    conversion was admitted after it.

    Every test here drives the REAL method. A constructor probe would have passed against
    the defect, because the constructor was never where the value was lost.
    """

    FAMILY = "Fake_Quarantined_Family"

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    def test_the_real_loader_attaches_the_entry_to_the_template_shape(self):
        """The half no test covered: the shape where the HIGH lived."""
        library = _auto_generated_library("a_seed", self.FAMILY)
        reaction = library.get_library_reactions()[0]

        assert isinstance(reaction, TemplateReaction)
        assert reaction.family == self.FAMILY
        assert getattr(reaction, "entry", None) is not None, (
            "the auto-generated template shape reached the model with no entry, so the "
            "longDesc that carries its authorship was unreachable")
        assert reaction.entry.long_desc.count("family:") == 1

    def test_a_seed_conversion_leaves_the_rate_refused(self, monkeypatch, tmp_path):
        """
        The acceptance for the HIGH, through `add_seed_mechanism_to_core`.

        The family is quarantined on disk and deliberately NOT loaded, which is exactly
        when the conversion fires.
        """
        _quarantined_on_disk(tmp_path)
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        library = _auto_generated_library("a_seed", QUARANTINED_ON_DISK)
        _model_database(monkeypatch, {"a_seed": library}, {})

        model = CoreEdgeReactionModel()
        with pytest.raises(QuarantinedKineticsError) as raised:
            model.add_seed_mechanism_to_core("a_seed")
        assert QUARANTINED_ON_DISK in str(raised.value)
        assert not model.core.reactions

    def test_a_reaction_library_conversion_leaves_the_rate_refused(self, monkeypatch,
                                                                   tmp_path):
        """The same defect at the second conversion site, through the real method."""
        _quarantined_on_disk(tmp_path)
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        library = _auto_generated_library("a_library", QUARANTINED_ON_DISK)
        _model_database(monkeypatch, {"a_library": library}, {})

        model = CoreEdgeReactionModel()
        with pytest.raises(QuarantinedKineticsError):
            model.add_reaction_library_to_edge("a_library")
        assert not model.edge.reactions

    def test_the_conversion_still_writes_the_library_label_into_the_family_slot(
            self, monkeypatch, tmp_path):
        """
        The round-89 addendum's contract, re-asserted at the site this round changed.

        `electron_placement.py` keys `FAMILY_ELECTRON_PLACEMENT` on `Reaction.family` and
        DEPENDS on a `LibraryReaction` putting its library label there. Carrying the entry
        across the conversion must not disturb that, so this drives a conversion that is
        allowed to complete and looks at the slot.
        """
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        library = _auto_generated_library("a_seed", "A_Family_Nobody_Quarantined")
        _model_database(monkeypatch, {"a_seed": library}, {})

        model = CoreEdgeReactionModel()
        model.add_seed_mechanism_to_core("a_seed")

        assert len(model.core.reactions) == 1
        converted = model.core.reactions[0]
        assert isinstance(converted, LibraryReaction)
        assert converted.family == "a_seed", (
            "the conversion stopped writing the library label into .family, which is the "
            "value electron_placement.py resolves argon ionisation through")
        assert authoring_families(converted) == ["A_Family_Nobody_Quarantined"], (
            "the authorship did not survive the conversion")


class TestTheDiskAnswerIsRevalidated:
    """
    Round 92. The disk cache stored one answer per (database, label) forever, with no
    check that the database still said it. Two consequences, both fail-open in one
    direction or the other, and one race that fails open specifically.
    """

    LABEL = "A_Family_That_Changes"

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    def _database(self, tmp_path, manifest=False):
        family = tmp_path / "kinetics" / "families" / self.LABEL
        family.mkdir(parents=True)
        if manifest:
            write_manifest(family)
        return family

    def test_a_manifest_added_after_the_first_lookup_is_seen(self, monkeypatch, tmp_path):
        family = self._database(tmp_path)
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))

        assert resolve_quarantine(self.LABEL) == (None, True)
        write_manifest(family)
        quarantine, answered = resolve_quarantine(self.LABEL)
        assert answered and quarantine is not None, (
            "the cache answered from a database state that no longer exists")

    def test_a_manifest_removed_after_the_first_lookup_is_seen(self, monkeypatch,
                                                              tmp_path):
        """The same staleness in the direction that keeps refusing a released family."""
        family = self._database(tmp_path, manifest=True)
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))

        assert resolve_quarantine(self.LABEL)[0] is not None
        os.remove(os.path.join(str(family), QUARANTINE_FILENAME))
        assert resolve_quarantine(self.LABEL) == (None, True)

    def test_a_manifest_that_vanishes_mid_read_is_unanswered_not_clean(
            self, monkeypatch, tmp_path):
        """
        The TOCTOU window, made deterministic.

        The file is removed from inside the check that decides the manifest is there,
        which is what a concurrent database edit does at an unpredictable moment. The
        wrong answer here is the dangerous one: "this family carries no manifest" is a
        clean bill of health, and the first version of the code cached it.

        The hook moved in round 95 and the reason is worth recording, because a race
        harness that no longer sits on the code path is a test that passes while
        measuring nothing. Round 92's window was between an `os.path.exists` and a later
        `open`; round 95 replaced that `open` with a descriptor-based read that never
        calls `exists`, and `_manifest_signature` with an `os.lstat`. Both stats are
        hooked, so the same test exercises the same window on either engine and its red
        state on an older one is a real answer rather than a harness mismatch.
        """
        from rmgpy.data.kinetics import quarantine as module

        family = self._database(tmp_path, manifest=True)
        manifest_path = os.path.join(str(family), QUARANTINE_FILENAME)
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))

        state = {"fired": False}

        def racing(real):
            def hook(path, *args, **kwargs):
                answer = real(path, *args, **kwargs)
                if str(path) == manifest_path and not state["fired"]:
                    state["fired"] = True
                    os.remove(manifest_path)      # the race, made deterministic
                return answer
            return hook

        monkeypatch.setattr(module.os, "lstat", racing(module.os.lstat))
        monkeypatch.setattr(module.os, "stat", racing(module.os.stat))
        assert state["fired"] is False
        answer = resolve_quarantine(self.LABEL)
        assert state["fired"] is True, (
            "the manifest was never removed, so this test did not exercise the race it "
            "is named for")
        assert answer == (None, False), (
            "a manifest that vanished mid-read was reported as a family with no manifest")
        assert not module._DISK_QUARANTINE_CACHE, (
            "the race was cached, so every later lookup inherits it")


class TestAFamilyLabelIsNotAPath:
    """
    Round 92. A family label arrives from a `family:` line in an entry's longDesc --
    user-authored text -- and was joined onto the database path unexamined. The file at
    the end of that path is then EXECUTED.
    """

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    def _planted(self, tmp_path):
        outside = tmp_path / "outside"
        outside.mkdir()
        write_manifest(outside)
        database = tmp_path / "db"
        (database / "kinetics" / "families" / "An_Ordinary_Family").mkdir(parents=True)
        return outside, database

    def test_an_absolute_label_is_refused(self, monkeypatch, tmp_path):
        outside, database = self._planted(tmp_path)
        monkeypatch.setitem(settings, "database.directory", str(database))

        assert resolve_quarantine(str(outside)) == (None, False), (
            "an absolute label made os.path.join discard the database prefix, and the "
            "manifest planted outside the database was executed")

    def test_a_dotdot_label_is_refused(self, monkeypatch, tmp_path):
        outside, database = self._planted(tmp_path)
        monkeypatch.setitem(settings, "database.directory", str(database))

        escape = os.path.join("..", "..", "..", "outside")
        assert resolve_quarantine(escape) == (None, False)

    def test_a_symlinked_family_pointing_outside_is_refused(self, monkeypatch, tmp_path):
        """
        What the syntactic check alone cannot catch: an ordinary-looking single-component
        label whose directory is a symlink out of the tree.
        """
        outside, database = self._planted(tmp_path)
        link = database / "kinetics" / "families" / "Looks_Ordinary"
        try:
            os.symlink(str(outside), str(link))
        except (OSError, NotImplementedError):  # pragma: no cover - platform dependent
            pytest.skip("this filesystem does not support symlinks")
        monkeypatch.setitem(settings, "database.directory", str(database))

        assert resolve_quarantine("Looks_Ordinary") == (None, False)

    def test_an_ordinary_label_still_resolves(self, monkeypatch, tmp_path):
        """The positive control: the refusal must not cost the ordinary case."""
        outside, database = self._planted(tmp_path)
        monkeypatch.setitem(settings, "database.directory", str(database))

        assert resolve_quarantine("An_Ordinary_Family") == (None, True)


CALLSITE_MANIFEST = """
name = "Synthetic/quarantine"
state = "QUARANTINED FOR TESTING"
appliesToKineticsClass = "Arrhenius"
reason = "a reason"
requiresEngineModule = "rmgpy.data.kinetics.quarantine"
requiresEngineSymbol = "check_quarantine"
requiresEngineCallSites = ("rmgpy.rmg.model",)
"""


def _module_carrying(tmp_path, source, tag):
    """
    A module object under the name the manifest declares, carrying `source`.

    Nothing is executed: the check reads the module's source and compares its binding of
    the symbol, so a module with `__file__` set and the real gate bound is exactly the
    shape it inspects. That is what lets these tests ask what the check accepts without
    editing `rmgpy/rmg/model.py` in the worktree.
    """
    import linecache
    import types

    from rmgpy.data.kinetics.quarantine import check_quarantine as real_gate

    path = tmp_path / "carrier_{0}.py".format(tag)
    path.write_text(source)
    linecache.checkcache(str(path))
    module = types.ModuleType("rmgpy.rmg.model")
    module.__file__ = str(path)
    module.check_quarantine = real_gate
    return module


def _load_with_module(tmp_path, tag, module):
    directory = tmp_path / ("manifest_" + tag)
    directory.mkdir()
    (directory / QUARANTINE_FILENAME).write_text(CALLSITE_MANIFEST)
    import sys

    previous = sys.modules.get("rmgpy.rmg.model")
    sys.modules["rmgpy.rmg.model"] = module
    try:
        return load_family_quarantine("Synthetic", str(directory))
    finally:
        if previous is not None:
            sys.modules["rmgpy.rmg.model"] = previous
        else:
            del sys.modules["rmgpy.rmg.model"]


class TestTheCallSitePinCoversEveryGate:
    """
    Round 92, and the third round this pin has been found weaker than its name.

    Round 87 found that `math.pi` satisfied the symbol check. Round 89 found that an
    import satisfied the call-site check. This round: a call the interpreter can never
    reach, a same-named call on an unrelated object, and -- the one that matters -- one
    surviving gate out of the four admission paths the engine actually has.

    The bound is still stated rather than claimed away: none of this proves a call
    EXECUTES. What it now proves is coverage of the named sites.
    """

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    DEAD = ("from rmgpy.data.kinetics.quarantine import check_quarantine\n"
            "\n"
            "\n"
            "class CoreEdgeReactionModel(object):\n"
            "    def add_reaction_to_core(self, rxn):\n"
            "        if False:\n"
            "            check_quarantine(rxn, stage='never runs')\n")

    UNRELATED = ("from rmgpy.data.kinetics.quarantine import check_quarantine\n"
                 "\n"
                 "\n"
                 "class CoreEdgeReactionModel(object):\n"
                 "    def add_reaction_to_core(self, rxn):\n"
                 "        somebody_else.check_quarantine(rxn)\n")

    def test_a_call_in_a_dead_branch_does_not_satisfy_it(self, tmp_path):
        with pytest.raises(DatabaseError):
            _load_with_module(tmp_path, "dead",
                              _module_carrying(tmp_path, self.DEAD, "dead"))

    def test_a_same_named_call_on_something_else_does_not_satisfy_it(self, tmp_path):
        with pytest.raises(DatabaseError):
            _load_with_module(tmp_path, "unrelated",
                              _module_carrying(tmp_path, self.UNRELATED, "unrelated"))

    def test_deleting_three_of_the_four_gates_is_refused(self, tmp_path):
        """
        The one a module-level call count cannot see. Three admission paths lose their
        gate, one keeps it, and the old check counted 1 and was satisfied.
        """
        import rmgpy.rmg.model as real_model

        source = inspect.getsource(real_model)
        thinned = source
        removed = 0
        for stage in ("admission to the model core",
                      "admission to the model edge",
                      "admission to a pressure-dependent network"):
            for line in source.split("\n"):
                if "check_quarantine(" in line and stage in line:
                    thinned = thinned.replace(line + "\n", "")
                    removed += 1
        assert removed == 3, "the gate calls this test deletes have moved"

        with pytest.raises(DatabaseError) as raised:
            _load_with_module(tmp_path, "thinned",
                              _module_carrying(tmp_path, thinned, "thinned"))
        message = str(raised.value)
        assert "add_reaction_to_core" in message
        assert "add_reaction_to_unimolecular_networks" in message

    def test_the_real_module_satisfies_it(self, tmp_path):
        """
        The positive control, and the tripwire for this engine growing a fifth admission
        path: the required list lives beside the checker, so a new ungated path fails
        here rather than in a run.
        """
        import rmgpy.rmg.model as real_model

        assert _load_with_module(tmp_path, "real", real_model) is not None


class TestTheSuppressionConsultsDiskToo:
    """
    Round 92. The `answered=False` policy admits a rate whose family this database does
    not have, and is defensible only because it is LOUD. The suppression that keeps that
    warning out of ordinary runs iterated LOADED families, so a run early enough that
    nothing is loaded yet went quiet in a database that does carry manifests.

    A correction to the review that raised it: an unloaded family which HAS a manifest on
    disk is refused, not admitted -- `resolve_quarantine` reads it. The silence is in the
    other case, a family the database does not contain at all, and it is the database's
    other manifests that should have kept the warning switched on.
    """

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    def _reaction(self):
        return make_library_reaction(
            library="a_foreign_seed",
            comment="family: A_Family_This_Database_Does_Not_Have")

    def test_a_database_whose_manifests_are_only_on_disk_still_warns(
            self, monkeypatch, tmp_path, caplog):
        families = tmp_path / "kinetics" / "families" / "A_Quarantined_Family"
        families.mkdir(parents=True)
        write_manifest(families)
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        _register_families(monkeypatch, {})

        reaction = self._reaction()
        with caplog.at_level(logging.WARNING):
            check_quarantine(reaction, stage="a test", kinetics=reaction.kinetics)

        assert any("Cannot tell whether" in record.getMessage()
                   for record in caplog.records), (
            "the rate was admitted in silence by a database that does quarantine "
            "something; the policy allowing the admission rests on the warning")

    def test_a_database_with_no_quarantine_anywhere_stays_quiet(
            self, monkeypatch, tmp_path, caplog):
        """
        The other half of the same policy, and the reason the suppression exists: an
        ordinary database loading a foreign seed must not get a line per reaction.
        """
        (tmp_path / "kinetics" / "families" / "An_Ordinary_Family").mkdir(parents=True)
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        _register_families(monkeypatch, {})

        reaction = self._reaction()
        with caplog.at_level(logging.WARNING):
            check_quarantine(reaction, stage="a test", kinetics=reaction.kinetics)

        assert not [record for record in caplog.records
                    if "Cannot tell whether" in record.getMessage()]


# ---------------------------------------------------------------------------------------
# Round 95
# ---------------------------------------------------------------------------------------

#: Quarantined on disk, deliberately not loaded.
QUARANTINED_FIRST = "A_Family_Quarantined_On_Disk"

#: Loaded, ordinary, unquarantined -- and named in a SECOND `family:` line. That is the
#: whole attack: the loader keeps the last label, and a loaded family is never converted
#: to a library reaction, so the shape round 89 fixed is never entered.
INNOCENT_SECOND = "An_Ordinary_Loaded_Family"


def _library_declaring(label, family_labels, electrons=0, degeneracy=1, auto=True,
                       elementary_high_p=False, allow_pdep_route=False,
                       allow_max_rate_violation=False):
    """
    An auto-generated library whose one entry declares `family_labels`, in order.

    `_auto_generated_library` declares exactly one. This one takes a list because the
    defect is what happens when a longDesc carries more than one, and carries `electrons`
    and `degeneracy` because the conversion used to drop both.

    The three flags go onto ``entry.item`` -- the blank `Reaction` -- which is exactly
    where `KineticsLibrary.load_entry` puts them when a library file declares them, and
    is the only place they can come from. Nothing here touches the loader's OUTPUT: the
    point of round 102's finding is that a test assigning these onto the built reaction
    asks the converter only to carry values the fixture put there.
    """
    long_desc = ["Matched reaction 3 Lip + CH3 <=> CH3Li in {0}/rate rule [Root]".format(
        family_labels[0]), "Euclidian distance = 0"]
    long_desc += ["family: {0}".format(name) for name in family_labels]
    library = KineticsLibrary(label=label, name=label)
    library.auto_generated = auto
    library.entries = {
        1: Entry(
            index=1,
            label="Lip + CH3 <=> CH3Li",
            item=Reaction(
                reactants=[Species(label="Lip", molecule=[Molecule(smiles="[Li+]")],
                                   reactive=False),
                           Species(label="CH3", molecule=[Molecule(smiles="[CH3]")],
                                   reactive=False)],
                products=[Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")],
                                  reactive=False)],
                reversible=False, electrons=electrons, degeneracy=degeneracy,
                elementary_high_p=elementary_high_p,
                allow_pdep_route=allow_pdep_route,
                allow_max_rate_violation=allow_max_rate_violation),
            data=_marcus_with(""),
            long_desc="\n".join(long_desc),
        )
    }
    return library


class TestASecondFamilyLineCannotHideTheFirstOnEitherShape:
    """
    Round 95's HIGH. `authoring_families`' docstring has promised since round 87 that a
    second ``family:`` line cannot hide the first, and that promise held for the library
    shape only: a template reaction took a short path and returned the single
    ``reaction.family`` slot, which the loader overwrites once per line.

    Round 92 attached the entry to this shape and these tests are what make the attachment
    load-bearing rather than decorative. The `add_reaction_to_core` test is the acceptance
    the manager named: the real loader, the real admission method, no constructor probe.
    """

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    def _on_disk(self, tmp_path):
        _quarantined_on_disk(tmp_path, QUARANTINED_FIRST)
        (tmp_path / "kinetics" / "families" / INNOCENT_SECOND).mkdir(parents=True)

    def test_authoring_families_reads_every_label_on_the_template_shape(
            self, monkeypatch, tmp_path):
        """The reader, on its own. Fixing only the loader would leave this open."""
        self._on_disk(tmp_path)
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        library = _library_declaring("a_seed", [QUARANTINED_FIRST, INNOCENT_SECOND])
        reaction = library.get_library_reactions()[0]

        assert isinstance(reaction, TemplateReaction)
        assert QUARANTINED_FIRST in authoring_families(reaction), (
            "the entry declares two families and the reader returned only the slot, so "
            "a quarantined first label is invisible behind an innocent second one")

    def test_the_hidden_label_is_refused_at_add_reaction_to_core(self, monkeypatch,
                                                                 tmp_path):
        """The acceptance: admitted at 9b89a937c, refused here, through the real method."""
        self._on_disk(tmp_path)
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        library = _library_declaring("a_seed", [QUARANTINED_FIRST, INNOCENT_SECOND])
        _model_database(monkeypatch, {"a_seed": library},
                        {INNOCENT_SECOND: _loaded_family(INNOCENT_SECOND)})
        reaction = library.get_library_reactions()[0]

        model = CoreEdgeReactionModel()
        with pytest.raises(QuarantinedKineticsError) as raised:
            model.add_reaction_to_core(reaction)
        assert QUARANTINED_FIRST in str(raised.value)
        assert not model.core.reactions

    def test_the_family_slot_is_not_normalised(self, monkeypatch, tmp_path):
        """
        The round-89 addendum, at the site this round changed. The repair is in the
        READER; the slot keeps holding the last label, because `model.py` compares it
        against the loaded families and `electron_placement.py` reads it.
        """
        self._on_disk(tmp_path)
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        library = _library_declaring("a_seed", [QUARANTINED_FIRST, INNOCENT_SECOND])
        reaction = library.get_library_reactions()[0]

        assert reaction.family == INNOCENT_SECOND, (
            "the slot's meaning was changed instead of the reader's, which the round-89 "
            "addendum forbids: the slot is what model.py and electron_placement.py read")

    def test_the_library_shape_still_refuses_both_labels(self, monkeypatch, tmp_path):
        """The guarantee that already held. It must keep holding."""
        self._on_disk(tmp_path)
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        library = _library_declaring("an_ordinary_library",
                                     [QUARANTINED_FIRST, INNOCENT_SECOND], auto=False)
        _model_database(monkeypatch, {"an_ordinary_library": library},
                        {INNOCENT_SECOND: _loaded_family(INNOCENT_SECOND)})
        reaction = library.get_library_reactions()[0]

        assert isinstance(reaction, LibraryReaction)
        with pytest.raises(QuarantinedKineticsError):
            model = CoreEdgeReactionModel()
            model.add_reaction_to_core(reaction)


def _loaded_family(label, quarantine=None, path=None):
    class _Family(object):
        pass

    family = _Family()
    family.label = label
    family.quarantine = quarantine
    family.quarantine_path = path
    return family


class TestTheConversionCarriesEveryField:
    """
    Round 95's MEDIUM 3. Both conversions constructed their replacement `LibraryReaction`
    with 8 of the 17 fields the class takes, so the rest reverted to defaults.

    `electrons` is the one with teeth in this campaign: a plasma reaction carries a signed
    electron count because that is what the charge-balance checks read, so
    ``Ar + e- => Ar+ + 2e-`` converted through this path claimed zero net electrons.
    """

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    def test_a_converted_electron_transfer_reaction_keeps_its_electrons(
            self, monkeypatch, tmp_path):
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        (tmp_path / "kinetics" / "families").mkdir(parents=True)
        library = _library_declaring("a_seed", ["A_Family_Nobody_Has"], electrons=1,
                                     degeneracy=3)
        _model_database(monkeypatch, {"a_seed": library}, {})

        model = CoreEdgeReactionModel()
        model.add_seed_mechanism_to_core("a_seed")

        assert len(model.core.reactions) == 1
        converted = model.core.reactions[0]
        assert isinstance(converted, LibraryReaction)
        assert converted.electrons == 1, (
            "the conversion reset the electron count to 0, so a seed reaction that "
            "produces a free electron arrived in the core claiming charge balance")
        assert converted.degeneracy == 3

    def test_every_field_the_class_takes_is_carried(self):
        """
        The enumeration, rather than a hand-written list of fields to check.

        The defect being closed IS a hand-written field list falling behind the class it
        builds, so a test that hard-codes its own list reproduces the defect one layer up.
        This one reads `LibraryReaction.__init__`'s signature, so a field added to the
        class later is covered without editing this file.
        """
        from rmgpy.rmg.model import as_library_reaction

        # Built by the real loader, from an entry that declares the flags the way a
        # library FILE declares them. Round 102's finding was that this test used to
        # assign `elementary_high_p`, `allow_pdep_route` and `allow_max_rate_violation`
        # onto the loader's output -- so the converter was only ever asked to carry
        # values the fixture had put there, and the loader dropping all three on the
        # template shape was invisible. Nothing below writes to `source`.
        source = _library_declaring(
            "a_seed", ["A_Family"], electrons=1, degeneracy=3,
            elementary_high_p=True, allow_pdep_route=True,
            allow_max_rate_violation=True).get_library_reactions()[0]

        assert isinstance(source, TemplateReaction)
        for flag in ("elementary_high_p", "allow_pdep_route", "allow_max_rate_violation"):
            assert getattr(source, flag) is True, (
                "the LOADER dropped {0!r} before the converter ever saw it; "
                "`elementary_high_p` lost here means the reaction silently misses "
                "pressure-dependent routing".format(flag))

        converted = as_library_reaction(source, "a_library")

        assert converted.library == "a_library"
        for name in inspect.signature(LibraryReaction.__init__).parameters:
            if name in ("self", "library", "reactants", "products"):
                continue
            was = getattr(source, name, None)
            assert getattr(converted, name, None) == was, (
                "the conversion dropped {0!r}: {1!r} became {2!r}".format(
                    name, was, getattr(converted, name, None)))


class TestTheManifestFileIsValidatedNotOnlyItsDirectory:
    """
    Round 95's path MEDIUM. Round 92 constrained the family *directory* to the database.
    The manifest inside it is what gets **executed**, and it could still be a symlink to
    anywhere; a label carrying a NUL raised an uncaught `ValueError` out of the gate, and
    one carrying a newline selected a real directory and executed its manifest.
    """

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    def _database(self, monkeypatch, tmp_path, label):
        family = tmp_path / "kinetics" / "families" / label
        family.mkdir(parents=True)
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        return family

    def test_a_manifest_symlinked_out_of_the_database_is_refused(self, monkeypatch,
                                                                 tmp_path):
        label = "A_Family_With_A_Symlinked_Manifest"
        family = self._database(monkeypatch, tmp_path, label)
        outside = tmp_path.parent / "not_a_manifest_{0}.py".format(tmp_path.name)
        outside.write_text(MANIFEST)
        os.symlink(str(outside), str(family / QUARANTINE_FILENAME))

        quarantine, answered = resolve_quarantine(label)
        assert quarantine is None and answered is False, (
            "the manifest was a link out of the database and it was followed and "
            "executed; refusing it must also leave the family UNanswered, since a "
            "refusal is not a clean bill of health")

    def test_a_label_carrying_a_nul_byte_is_refused_rather_than_raising(self, monkeypatch,
                                                                       tmp_path):
        self._database(monkeypatch, tmp_path, "A_Family")
        assert resolve_quarantine("A_Family\x00.py") == (None, False)

    def test_a_label_carrying_a_newline_is_refused(self, monkeypatch, tmp_path):
        label = "A_Family\nOhNo"
        try:
            family = self._database(monkeypatch, tmp_path, label)
        except OSError:
            pytest.skip("this filesystem will not hold a newline in a directory name")
        write_manifest(str(family))

        assert resolve_quarantine(label) == (None, False), (
            "a newline in a family label selected a real directory and its manifest "
            "was executed")

    def test_an_ordinary_label_still_resolves(self, monkeypatch, tmp_path):
        label = "An_Ordinary_Family"
        family = self._database(monkeypatch, tmp_path, label)
        write_manifest(str(family))

        quarantine, answered = resolve_quarantine(label)
        assert answered and quarantine is not None, (
            "the hardening refused an ordinary manifest; the check has stopped being a "
            "check and started being a wall")

    def test_a_manifest_that_is_not_a_regular_file_is_refused(self, monkeypatch,
                                                              tmp_path):
        """A directory named quarantine.py would otherwise raise from inside the gate."""
        label = "A_Family_Whose_Manifest_Is_A_Directory"
        family = self._database(monkeypatch, tmp_path, label)
        (family / QUARANTINE_FILENAME).mkdir()

        assert resolve_quarantine(label) == (None, False)


class TestALoadedFamilyRereadsItsManifest:
    """
    Round 95's cache MEDIUM. Round 92 made the *unloaded* half re-stat the manifest on
    every lookup. The loaded half returned `family.quarantine` before it looked at disk,
    so a manifest added, edited or removed after load was invisible for the life of the
    process -- the same staleness, on the other branch of the same function.
    """

    LABEL = "A_Loaded_Family_Whose_Manifest_Appears"

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    def _loaded(self, monkeypatch, tmp_path, quarantine=None):
        family_dir = tmp_path / "kinetics" / "families" / self.LABEL
        family_dir.mkdir(parents=True)
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        _register_families(monkeypatch, {self.LABEL: quarantine})
        return family_dir

    def test_a_manifest_added_after_the_family_loaded_is_seen(self, monkeypatch,
                                                              tmp_path):
        family_dir = self._loaded(monkeypatch, tmp_path)
        assert resolve_quarantine(self.LABEL) == (None, True)

        write_manifest(str(family_dir))
        quarantine, answered = resolve_quarantine(self.LABEL)
        assert answered and quarantine is not None, (
            "the family was loaded before the manifest existed and the gate went on "
            "answering from the object for the life of the process")

    def test_a_manifest_removed_after_the_family_loaded_is_seen(self, monkeypatch,
                                                                tmp_path):
        family_dir = self._loaded(monkeypatch, tmp_path)
        write_manifest(str(family_dir))
        assert resolve_quarantine(self.LABEL)[0] is not None

        os.remove(str(family_dir / QUARANTINE_FILENAME))
        assert resolve_quarantine(self.LABEL) == (None, True)


class TestTheSuppressionRescansRatherThanTrustingACachedNo:
    """
    Round 95's second cache MEDIUM. `_DISK_ANY_QUARANTINE_CACHE` is keyed on the families
    ROOT's signature, and writing a manifest inside an existing family directory does not
    move it -- so a cached "this database quarantines nothing" outlived the fact, and
    `_warn_unanswered` stayed silent in a database that had started quarantining something.
    """

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    def test_a_manifest_added_inside_an_existing_family_ends_the_suppression(
            self, monkeypatch, tmp_path):
        from rmgpy.data.kinetics.quarantine import _database_has_any_quarantine

        family = tmp_path / "kinetics" / "families" / "A_Family_Gaining_A_Manifest"
        family.mkdir(parents=True)
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        _register_families(monkeypatch, {})

        assert _database_has_any_quarantine() is False
        write_manifest(str(family))
        assert _database_has_any_quarantine() is True, (
            "the negative was cached against the families root's mtime, which a write "
            "inside an existing family does not move")


class TestThePinDocumentsWhatItActuallyChecks:
    """
    Round 95's fourth finding on `requiresEngineCallSites`, and the option taken.

    The manager's ruling was "either it verifies what its name says, or rename it and
    correct the docstring to state exactly what it checks". The check is renamed
    internally to `_gate_calls_in_source` and its documentation now enumerates what it
    does NOT verify. It is deliberately NOT strengthened a fifth time.
    """

    def test_the_documentation_enumerates_what_it_does_not_verify(self):
        from rmgpy.data.kinetics import quarantine as module

        checker = getattr(module, "_gate_calls_in_source", None)
        assert checker is not None, (
            "the check still carries a name that promises reachability analysis")
        text = ((checker.__doc__ or "")
                + (module._check_engine_requirements.__doc__ or "")).lower()
        for missing in ("execut", "dominat", "argument", "propagat", "short-circuit",
                        "shadow", "lambda"):
            assert missing in text, (
                "the documentation does not say that the check fails to verify {0!r}; "
                "a syntactic check whose prose implies more is the defect four rounds "
                "have found".format(missing))


class TestTheCarrierSurvivesATransform:
    """
    Round 95's LOW. `get_library_reactions` attaches the entry, and
    `TemplateReaction.__reduce__` and `copy()` enumerate their fields by hand -- neither
    listed it, so a pickle or a copy dropped the provenance. The same carrier-loss class
    as the HIGH, one transform over.
    """

    def _built(self):
        library = _library_declaring("a_seed", ["A_Family"])
        return library.get_library_reactions()[0]

    def test_the_entry_survives_a_pickle(self):
        reaction = self._built()
        assert reaction.entry is not None
        assert pickle.loads(pickle.dumps(reaction)).entry is not None

    def test_the_entry_survives_a_copy(self):
        reaction = self._built()
        assert reaction.entry is not None
        assert reaction.copy().entry is not None


class TestTheDirectoryCannotBeSwappedUnderTheCheck:
    """
    Round 99. Round 95 validated the family directory and then opened it *by name*, so the
    check and the use were two independent resolutions of the same string with a window
    between them. A directory symlink swapped into that window was followed, and the file
    at the end of it is `exec()`d.

    The closure is structural rather than another check: the label is opened as a single
    component through the families root's own descriptor with `O_NOFOLLOW`, and the
    manifest through *that* descriptor. Nothing is re-derived from a string after being
    checked, so there is no window left to race.

    The price is that a family directory which is itself a symlink is refused even when it
    points inside the database. Measured before taking it: zero symlinks anywhere under
    `input/kinetics/families` in either database on this box.
    """

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    def _database(self, monkeypatch, tmp_path, label):
        family = tmp_path / "kinetics" / "families" / label
        family.mkdir(parents=True)
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        return family

    def test_a_directory_swapped_during_the_check_is_not_followed(self, monkeypatch,
                                                                  tmp_path):
        """
        The window made deterministic, by swapping from inside the check itself.

        `_family_directory` calls `os.path.realpath` on the family path; that call IS the
        containment check. It is allowed to compute the true, in-tree answer, and the
        directory is replaced with a link out of the database before that answer is
        returned. The check therefore passes on the directory that was there, and anything
        that reaches for the *name* afterwards gets the one that is there now.
        """
        label = "A_Family_Swapped_Under_The_Check"
        family = self._database(monkeypatch, tmp_path, label)
        outside = tmp_path.parent / "outside_{0}".format(tmp_path.name)
        outside.mkdir()
        write_manifest(str(outside))

        state = {"fired": False}
        original = os.path.realpath
        target = os.path.abspath(str(family))

        def realpath(path, *args, **kwargs):
            answer = original(path, *args, **kwargs)
            if not state["fired"] and os.path.abspath(path) == target:
                state["fired"] = True
                os.rename(target, target + ".was_here")
                os.symlink(str(outside), target)
            return answer

        monkeypatch.setattr(os.path, "realpath", realpath)
        quarantine, answered = resolve_quarantine(label)

        assert state["fired"] is True, (
            "this test did not exercise the race it is named for: the containment check "
            "never resolved the family path, so nothing was swapped")
        assert quarantine is None and answered is False, (
            "a directory swapped between the containment check and the open was followed, "
            "and a manifest from outside the database was executed; the refusal must also "
            "leave the family UNanswered, since a refusal is not a clean bill of health")

    def test_a_family_directory_that_is_a_symlink_is_refused_even_inside_the_tree(
            self, monkeypatch, tmp_path):
        """
        The deliberate narrowing, pinned so it is not quietly relaxed back.

        Admitting an in-tree link means containment is decided by resolving a name a
        second time, which is the defect. A family directory is a real directory, a
        direct child of the families root -- that rule is checkable in one sentence and
        has no window in it.
        """
        real = "A_Real_Family"
        link = "A_Family_Linked_To_Its_Sibling"
        family = self._database(monkeypatch, tmp_path, real)
        write_manifest(str(family))
        families = tmp_path / "kinetics" / "families"
        os.symlink(str(family), str(families / link))

        assert resolve_quarantine(link) == (None, False), (
            "a family directory that is a symbolic link was resolved and its manifest "
            "executed; it is refused, and refused as unanswered")
        assert resolve_quarantine(real)[1] is True, (
            "the real family behind the link must still answer -- the refusal is of the "
            "link, not of what it points at")

    def test_a_database_root_that_is_itself_a_symlink_still_resolves(self, monkeypatch,
                                                                     tmp_path):
        """
        The narrowing stops at the label. Everything above it is configuration.

        `database.directory` is a local setting, not something a shipped `family:` line
        can steer, and symlinking a whole checkout is ordinary. Refusing links all the way
        up would turn this check into a wall for a case no entry can reach.
        """
        real_root = tmp_path / "real_root"
        label = "A_Family_Under_A_Symlinked_Root"
        family = real_root / "kinetics" / "families" / label
        family.mkdir(parents=True)
        write_manifest(str(family))
        linked_root = tmp_path / "linked_root"
        os.symlink(str(real_root), str(linked_root))
        monkeypatch.setitem(settings, "database.directory", str(linked_root))

        quarantine, answered = resolve_quarantine(label)
        assert answered and quarantine is not None, (
            "a database reached through a symlinked root stopped resolving its own "
            "manifests; the pin belongs on the label, not on the path to the database")


class TestAnUnreadableManifestIsNotAnAbsentOne:
    """
    Round 99's HIGH. `_manifest_signature` mapped **every** `OSError` from its `lstat` to
    `'absent'` -- permission denied, an I/O error, a descriptor exhaustion -- and
    `resolve_quarantine` turns "the family directory is there and carries no manifest"
    into `(None, True)`: a clean bill of health, then cached for the life of the run.

    This is round 89's HIGH 2 one layer down. That round separated "authorship
    unrecoverable" from "manifest never consulted"; this path collapsed "there is no
    manifest" into "I could not look". A manifest that exists and cannot be read is an
    unanswered question and takes the unanswered path, uncached.

    The demonstration is a permission error rather than a deleted file, because a deleted
    file genuinely is absent and would prove nothing about the collapse.
    """

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    def _unreadable_family(self, monkeypatch, tmp_path, label):
        family = tmp_path / "kinetics" / "families" / label
        family.mkdir(parents=True)
        write_manifest(str(family))
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        return _deny_access(monkeypatch, family)

    def test_a_manifest_that_cannot_be_examined_is_unanswered_not_clean(self, monkeypatch,
                                                                        tmp_path):
        label = "A_Family_Whose_Manifest_Cannot_Be_Read"
        self._unreadable_family(monkeypatch, tmp_path, label)

        assert resolve_quarantine(label) == (None, False), (
            "a manifest that exists and could not be examined was reported as a family "
            "carrying no manifest at all -- a clean bill of health handed out by the "
            "check that failed")

    def test_the_unreadable_answer_is_not_entered_into_the_cache(self, monkeypatch,
                                                                 tmp_path):
        """
        An answer produced by a failure is not a property of the cache key.

        Measured limit, stated because the round's brief claimed more: the wrong answer
        did **not** in fact survive the condition clearing, because restoring the
        permission moves the signature from `(True, None, 'absent')` to
        `(True, <identity>, 'regular')` and misses the cache anyway. What is pinned here
        is the narrower true thing -- the key is not written at all -- which is the same
        rule round 95 arrived at for the loaded-family fallback.
        """
        from rmgpy.data.kinetics import quarantine as module

        label = "A_Family_Readable_Again_Later"
        allow = self._unreadable_family(monkeypatch, tmp_path, label)

        assert resolve_quarantine(label) == (None, False)
        assert not [key for key in module._DISK_QUARANTINE_CACHE if key[1] == label], (
            "the answer taken from a permission error was cached, under a key naming a "
            "directory and a label -- neither of which is what produced it")

        allow()
        second = resolve_quarantine(label)
        assert second[1] is True and second[0] is not None, (
            "the manifest stayed invisible after it became readable")


class TestTheSignatureNoticesAnEqualLengthEdit:
    """
    Round 99's cache-identity MEDIUM. The identity was `(mtime_ns, size, inode)`, and an
    in-place edit of the same length with the mtime put back preserves all three -- so a
    manifest whose `state` or `reason` changed kept answering with the old text for the
    rest of the run.

    `st_ctime_ns` closes it: the inode change time moves on any write and, unlike mtime,
    no userspace call can set it back. It costs nothing -- the `lstat` is already taken.
    Hashing the content would be stronger still and is deliberately not done: it turns
    every lookup of every family into a file read, for a case ctime already covers
    everywhere ctime is real.
    """

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    def test_an_in_place_edit_of_the_same_length_with_mtime_restored_is_seen(
            self, monkeypatch, tmp_path):
        label = "A_Family_Edited_In_Place"
        family = tmp_path / "kinetics" / "families" / label
        family.mkdir(parents=True)
        manifest = family / QUARANTINE_FILENAME
        before = MANIFEST.replace("a reason that must reach the error message",
                                  "reason A -- xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        after = MANIFEST.replace("a reason that must reach the error message",
                                 "reason B -- xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        assert len(before) == len(after), "the edit must not change the file's length"
        manifest.write_text(before)
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))

        first = resolve_quarantine(label)
        assert first[0] is not None and first[0].reason.startswith("reason A")

        stamps = os.stat(str(manifest))
        with open(str(manifest), "r+") as handle:
            handle.seek(0)
            handle.write(after)
        os.utime(str(manifest), ns=(stamps.st_atime_ns, stamps.st_mtime_ns))
        repeat = os.stat(str(manifest))
        assert (repeat.st_mtime_ns, repeat.st_size, repeat.st_ino) == (
            stamps.st_mtime_ns, stamps.st_size, stamps.st_ino), (
            "this test did not exercise the case it is named for: the edit moved one of "
            "the three fields the old identity was built from")

        second = resolve_quarantine(label)
        assert second[0] is not None and second[0].reason.startswith("reason B"), (
            "an edited manifest kept answering with its old text: the cache identity is "
            "preserved by an equal-length in-place edit with the mtime put back")


class TestTheAffectedSetCannotPassVacuously:
    """
    Round 99's third test-quality MEDIUM. `affected_entries` returns two lists, and every
    assertion made about their contents is an `all(...)` -- which is true of nothing.
    Neither branch of the enumeration had ever been driven to empty, so nothing showed
    what the guard against an empty result is.
    """

    def _family(self, rules, training):
        class _Depository:
            def __init__(self, label, entries):
                self.label = label
                self.entries = entries

        class _Rules:
            def __init__(self, entries):
                self.entries = entries

        class _Family:
            def __init__(self):
                self.rules = _Rules(rules)
                self.depositories = [_Depository("X/training", training)]

        return _Family()

    def _entry(self, data):
        return Entry(index=1, label="e", data=data)

    def test_a_family_with_no_matching_entries_enumerates_nothing(self, quarantine):
        affected = quarantine.affected_entries(self._family({}, {}))
        assert affected == {"rules": [], "training": []}
        assert all(isinstance(e.data, Marcus)
                   for e in affected["rules"] + affected["training"]), (
            "this is the vacuous pass itself, asserted so it is on the record: `all()` "
            "over the empty enumeration is true, which is why a count or a non-emptiness "
            "check has to sit beside it")

    def test_the_enumeration_is_selective_within_each_of_the_two_branches(self, quarantine):
        matching = self._entry(make_marcus())
        other = self._entry(Arrhenius(A=(1.0, "cm^3/(mol*s)"), n=0, Ea=(0, "kJ/mol")))
        family = self._family({"n": [matching, other]}, {"a": matching, "b": other})
        affected = quarantine.affected_entries(family)
        assert len(affected["rules"]) == 1 and affected["rules"][0] is matching
        assert len(affected["training"]) == 1 and affected["training"][0] is matching

    def test_a_family_with_rules_but_no_training_depository_is_not_an_error(self, quarantine):
        class _Family:
            rules = None
            depositories = []

        affected = quarantine.affected_entries(_Family())
        assert affected == {"rules": [], "training": []}


class TestTheCallSiteFieldIsNamedForWhatItChecks:
    """
    Round 99's rename. Five rounds have each produced one more construct that satisfies
    `requiresEngineCallSites` while gating nothing -- `math.pi`, import-without-call, a
    dead `if False:`, `False and check_quarantine(r)`, and a locally shadowed name. The
    ruling was not to strengthen it a sixth time but to make the field say what it does:
    a static syntactic presence check.

    The old spelling stays readable on purpose. Manifests declaring it ship in database
    repositories this engine change may not edit, and a renamed field that silently stops
    being read disarms a pin nobody is watching -- which is the failure this whole field
    group exists to prevent. Declaring both is refused rather than resolved by guessing.
    """

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    def _manifest(self, tmp_path, field):
        return write_manifest(
            str(tmp_path),
            MANIFEST
            + 'requiresEngineModule = "rmgpy.data.kinetics.quarantine"\n'
            + 'requiresEngineSymbol = "check_quarantine"\n'
            + '{0} = ("rmgpy.rmg.model",)\n'.format(field))

    def test_the_new_spelling_is_honoured(self, tmp_path):
        path = self._manifest(tmp_path, "requiresEngineCallsInSource")
        assert load_family_quarantine("Fake_Quarantined_Family", path) is not None

    def test_the_new_spelling_still_refuses_a_module_that_never_calls_the_gate(
            self, tmp_path, monkeypatch):
        """The rename must not have turned the check into a no-op under its new name."""
        module_dir = tmp_path / "fake_site"
        module_dir.mkdir()
        (module_dir / "binds_but_never_calls.py").write_text(
            "from rmgpy.data.kinetics.quarantine import check_quarantine\n")
        monkeypatch.syspath_prepend(str(module_dir))

        manifest_dir = tmp_path / "manifest"
        manifest_dir.mkdir()
        write_manifest(str(manifest_dir),
                       MANIFEST
                       + 'requiresEngineModule = "rmgpy.data.kinetics.quarantine"\n'
                       + 'requiresEngineSymbol = "check_quarantine"\n'
                       + 'requiresEngineCallsInSource = ("binds_but_never_calls",)\n')
        with pytest.raises(DatabaseError) as exc:
            load_family_quarantine("Fake_Quarantined_Family", str(manifest_dir))
        assert "never calls it" in str(exc.value)

    def test_the_old_spelling_is_still_honoured_and_says_so(self, tmp_path, caplog):
        path = self._manifest(tmp_path, "requiresEngineCallSites")
        with caplog.at_level(logging.WARNING):
            assert load_family_quarantine("Fake_Quarantined_Family", path) is not None
        assert "requiresEngineCallsInSource" in caplog.text, (
            "a manifest using the old spelling was read without being told the new one; "
            "a rename nobody is informed of is a rename that never happens in the data")

    def test_the_old_spelling_still_refuses_what_it_always_refused(self, tmp_path):
        """The alias must carry the check, not merely be accepted."""
        body = (MANIFEST
                + 'requiresEngineModule = "rmgpy.data.kinetics.quarantine"\n'
                + 'requiresEngineSymbol = "check_quarantine"\n'
                + 'requiresEngineCallSites = ("os",)\n')
        with pytest.raises(DatabaseError) as exc:
            load_family_quarantine("Fake_Quarantined_Family",
                                   write_manifest(str(tmp_path), body))
        assert "wired into" in str(exc.value)

    def test_declaring_both_spellings_is_refused(self, tmp_path):
        body = (MANIFEST
                + 'requiresEngineModule = "rmgpy.data.kinetics.quarantine"\n'
                + 'requiresEngineSymbol = "check_quarantine"\n'
                + 'requiresEngineCallsInSource = ("rmgpy.rmg.model",)\n'
                + 'requiresEngineCallSites = ("rmgpy.rmg.model",)\n')
        with pytest.raises(DatabaseError) as exc:
            load_family_quarantine("Fake_Quarantined_Family",
                                   write_manifest(str(tmp_path), body))
        assert "declare one" in str(exc.value)

    def test_the_documentation_of_the_field_states_what_it_actually_checks(self):
        """
        Asserted positively, on the words that must be PRESENT.

        Round 95 pinned the sibling docstring by checking that a false phrase was absent,
        and the corrected text contained that phrase inside the sentence disowning it --
        an assertion satisfied for the wrong reason. The content wanted here is a name
        that matches the check and a sentence that states its limit, so that is what is
        asserted.
        """
        from rmgpy.data.kinetics.quarantine import _check_engine_requirements

        text = _check_engine_requirements.__doc__
        assert "requiresEngineCallsInSource" in text
        assert "static" in text and "syntactic" in text, (
            "the field's own documentation must say it is a static syntactic check; "
            "four rounds of overstatement came from prose that implied more")
        for absent in ("executes", "dominates", "arguments", "propagates"):
            assert absent in text, (
                "the documentation must enumerate what the check does NOT verify -- it "
                "does not mention {0!r}".format(absent))


class TestEveryComponentBelowTheAnchorIsPinned:
    """
    Round 101's HIGH 1. Round 99 pinned the label component and the manifest, and opened
    the *parent* -- `<database>/kinetics/families` -- following symlinks. So the descent
    protected its last two components and not the path.

    `_family_directory`'s containment check cannot catch this, and that is the instructive
    part: it compares `realpath(family_path)` against `realpath(families_root)`, and when
    `families` is itself a link BOTH resolve outside the database, so the prefix test
    passes and reports containment in a tree that is no longer the database.

    There is always exactly one anchor that must be followed, because it may legitimately
    be a link -- here it is `settings['database.directory']`, which is local configuration.
    Everything below it is opened with `O_NOFOLLOW`. That is the whole rule.
    """

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    def _escape_at(self, monkeypatch, tmp_path, component):
        """Replace `component` of the descent with a link to a family outside the tree."""
        label = "A_Family_Reached_Through_A_Linked_Parent"
        outside = tmp_path / "outside_the_database"
        (outside / label).mkdir(parents=True)
        write_manifest(str(outside / label))

        root = tmp_path / "database"
        if component == "families":
            (root / "kinetics").mkdir(parents=True)
            os.symlink(str(outside), str(root / "kinetics" / "families"))
        else:
            root.mkdir()
            os.symlink(str(outside), str(root / "kinetics"))
            # `kinetics` links to a directory holding the family directly, so the
            # families level has to exist inside it for the join to land.
            (outside / "families").mkdir()
            os.rename(str(outside / label), str(outside / "families" / label))
        monkeypatch.setitem(settings, "database.directory", str(root))
        return label

    def test_a_symlinked_families_directory_does_not_reach_outside(self, monkeypatch,
                                                                   tmp_path):
        label = self._escape_at(monkeypatch, tmp_path, "families")
        assert resolve_quarantine(label) == (None, False), (
            "`kinetics/families` was a link out of the database and the manifest behind "
            "it was read and executed; the containment check cannot see this, because "
            "with that link in place the root and the family resolve to the same foreign "
            "tree and the prefix test passes")

    def test_a_symlinked_kinetics_directory_does_not_reach_outside(self, monkeypatch,
                                                                    tmp_path):
        label = self._escape_at(monkeypatch, tmp_path, "kinetics")
        assert resolve_quarantine(label) == (None, False), (
            "the escape moves one level up and must be refused at every level, not at "
            "the one the last repair happened to name")

    def test_the_configured_database_root_may_still_be_a_link(self, monkeypatch, tmp_path):
        """
        The anchor, and why there has to be one.

        A path cannot be descended without something to descend *from*, and that first
        thing is always followed. Choosing `database.directory` puts it on local
        configuration rather than on anything a shipped `family:` line can steer.
        """
        label = "A_Family_Under_A_Linked_Root"
        real_root = tmp_path / "real_root"
        (real_root / "kinetics" / "families" / label).mkdir(parents=True)
        write_manifest(str(real_root / "kinetics" / "families" / label))
        os.symlink(str(real_root), str(tmp_path / "linked_root"))
        monkeypatch.setitem(settings, "database.directory", str(tmp_path / "linked_root"))

        quarantine, answered = resolve_quarantine(label)
        assert answered and quarantine is not None, (
            "pinning every component including the anchor would refuse an ordinary "
            "symlinked checkout, which no shipped entry can influence")


class TestContainmentIsRefusedRatherThanApproximated:
    """
    Round 101's HIGH 2. The descent had an `else` branch that opened the whole pathname
    when directory descriptors were unavailable **or the last component was empty**, and
    said so in a comment. A documented hole is still a hole, and the empty case was
    reachable from the public loader.
    """

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    def test_an_empty_family_path_does_not_read_the_working_directory(self, tmp_path,
                                                                      monkeypatch):
        """
        `os.path.split('')` gives `('', '')`, and `os.path.join('', 'quarantine.py')` is
        `'quarantine.py'` -- a relative name. The fallback opened it, so an empty path
        executed whatever `quarantine.py` happened to sit in the process's working
        directory.
        """
        write_manifest(str(tmp_path))
        monkeypatch.chdir(str(tmp_path))

        from rmgpy.data.kinetics.quarantine import QUARANTINE_UNREADABLE

        answer = load_family_quarantine("A_Family", "")
        assert not isinstance(answer, KineticsQuarantine), (
            "an empty family path read `quarantine.py` from the current working "
            "directory and executed it")
        assert answer is QUARANTINE_UNREADABLE, (
            "round 111: the refusal must not be spelled `None`, which is what every "
            "ordinary family's clean bill of health is spelled")

    def test_a_root_family_path_is_refused(self, tmp_path):
        from rmgpy.data.kinetics.quarantine import QUARANTINE_UNREADABLE

        answer = load_family_quarantine("A_Family", os.sep)
        assert not isinstance(answer, KineticsQuarantine)
        assert answer is QUARANTINE_UNREADABLE

    def test_the_descent_is_refused_when_directory_descriptors_are_unavailable(
            self, monkeypatch, tmp_path):
        """
        The platform requirement, enforced rather than documented.

        Without `dir_fd` there is no way to open a component relative to a descriptor, so
        containment cannot be guaranteed at all -- and the previous fallback answered
        anyway. Refusing means a hypothetical platform loses the quarantine gate and is
        told so, which is the failure mode this module exists to prefer.
        """
        label = self._escape_at(monkeypatch, tmp_path, "families")
        monkeypatch.setattr(os, "supports_dir_fd", set())

        assert resolve_quarantine(label) == (None, False), (
            "with no directory descriptors the code opened the joined pathname and "
            "followed every link in it, reaching a manifest outside the database; there "
            "is no containment available on that path, so it must refuse rather than "
            "answer")

    _escape_at = TestEveryComponentBelowTheAnchorIsPinned._escape_at


class TestAnUnreadableDirectoryIsUnansweredForALoadedFamilyToo:
    """
    Round 101's HIGH 3 -- round 99's HIGH surviving in the branch round 99 did not test.

    `_manifest_signature` reports `(False, None, 'unreadable')` when the family directory
    cannot be examined *and* `os.path.isdir` cannot see it either, which is what happens
    when the permission error is on `kinetics/families` rather than on the family itself.
    `resolve_quarantine` then tested `directory_exists` before `kind`, so a **loaded**
    family took the "not where this database would put it" branch and returned its own
    attribute with `answered=True`. For a family whose attribute is `None` that is a clean
    bill of health produced by a permission error.

    The round-99 test could not see it: it registered no loaded family, and it made the
    family directory unreadable rather than its parent.
    """

    LABEL = "A_Loaded_Family_Behind_An_Unreadable_Parent"

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    def _unreadable_parent(self, monkeypatch, tmp_path, quarantine):
        families = tmp_path / "kinetics" / "families"
        (families / self.LABEL).mkdir(parents=True)
        write_manifest(str(families / self.LABEL))
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        _register_families(monkeypatch, {self.LABEL: quarantine})
        _deny_access(monkeypatch, families)

    def test_a_loaded_family_with_no_quarantine_is_unanswered_not_clean(self, monkeypatch,
                                                                        tmp_path):
        self._unreadable_parent(monkeypatch, tmp_path, None)
        assert resolve_quarantine(self.LABEL) == (None, False), (
            "a permission error on the families directory was reported as a loaded "
            "family that carries no quarantine -- a clean bill of health issued by the "
            "failure of the check, on the branch the last repair did not cover")

    def test_a_loaded_family_that_does_carry_one_is_also_unanswered(self, monkeypatch,
                                                                     tmp_path, quarantine):
        """
        The object's answer is not wrong here, and it is still refused.

        Answering from the attribute is exactly the staleness round 95 removed: it is
        whatever was true at load time, and the question asked is what is true now. When
        the disk cannot be consulted the honest answer is that there is no answer.
        """
        self._unreadable_parent(monkeypatch, tmp_path, quarantine)
        assert resolve_quarantine(self.LABEL) == (None, False)


class TestDeclaringBothSpellingsIsRefusedEvenWhenOneIsEmpty:
    """
    Round 101's LOW. The ambiguity check tested truthiness, not presence, so
    `requiresEngineCallsInSource = ()` beside a populated `requiresEngineCallSites` was
    not "both declared" -- and the promise made in the refusal message was false for
    exactly the manifest most likely to be written during a rename.
    """

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    GATE = ('requiresEngineModule = "rmgpy.data.kinetics.quarantine"\n'
            'requiresEngineSymbol = "check_quarantine"\n')

    @pytest.mark.parametrize("new_value,old_value", [
        ("()", '("rmgpy.rmg.model",)'),
        ('("rmgpy.rmg.model",)', "()"),
        ("()", "()"),
    ])
    def test_both_fields_present_is_refused_whatever_they_contain(self, tmp_path,
                                                                  new_value, old_value):
        body = (MANIFEST + self.GATE
                + 'requiresEngineCallsInSource = {0}\n'.format(new_value)
                + 'requiresEngineCallSites = {0}\n'.format(old_value))
        with pytest.raises(DatabaseError) as exc:
            load_family_quarantine("Fake_Quarantined_Family",
                                   write_manifest(str(tmp_path), body))
        assert "declare one" in str(exc.value)

    def test_an_empty_declaration_of_the_new_field_alone_is_still_a_declaration(
            self, tmp_path):
        """
        Presence, not truthiness, on the single-field path too: declaring the field with
        no modules pins nothing, and it must not be read as "the field is absent".
        """
        body = (MANIFEST + self.GATE + 'requiresEngineCallsInSource = ()\n')
        with pytest.raises(DatabaseError) as exc:
            load_family_quarantine("Fake_Quarantined_Family",
                                   write_manifest(str(tmp_path), body))
        assert "names no module" in str(exc.value)


class TestTheSignatureAndTheContentComeFromOneOpen:
    """
    Round 102's HIGH 1 -- round 99's lesson one level up.

    `_manifest_signature` stats the path; `load_family_quarantine` then opens the same
    path independently. A rename between those two resolutions reads B and stores B's
    quarantine **under A's signature**, after which restoring A is a cache HIT that
    returns B's answer -- A's criterion bypassed without A ever being read. `st_ctime_ns`
    cannot help, because the cached one is A's.

    The repair is not another check. The identity is taken with `os.fstat` on the
    descriptor the bytes came from, so the key and the value describe one object.
    """

    LABEL = "A_Family_Whose_Manifest_Is_Renamed_Mid_Lookup"

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    def _bodies(self):
        return (MANIFEST.replace("a reason that must reach the error message", "manifest A"),
                MANIFEST.replace("a reason that must reach the error message", "manifest B"))

    def test_a_rename_between_the_stat_and_the_read_does_not_outlive_itself(
            self, monkeypatch, tmp_path):
        family = tmp_path / "kinetics" / "families" / self.LABEL
        family.mkdir(parents=True)
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        a, b = self._bodies()
        manifest = family / QUARANTINE_FILENAME
        manifest.write_text(a)
        other = family / "manifest_b.py"
        other.write_text(b)

        state = {"fired": False}
        original = os.lstat
        target = os.path.abspath(str(manifest))

        def lstat(path, *args, **kwargs):
            answer = original(path, *args, **kwargs)
            if not state["fired"] and os.path.abspath(str(path)) == target:
                state["fired"] = True
                os.rename(target, str(family / "manifest_a.py"))
                os.rename(str(other), target)
            return answer

        from rmgpy.data.kinetics import quarantine as module

        monkeypatch.setattr(os, "lstat", lstat)
        first = resolve_quarantine(self.LABEL)
        monkeypatch.undo()
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))

        assert state["fired"] is True, (
            "this test did not exercise the rename it is named for: the signature was "
            "never taken on the manifest path")
        assert first[0].reason == "manifest B", (
            "the read is supposed to have picked up B -- that is the race, not the "
            "defect; the defect is what the cache now believes")

        # `quarantine.py` is B now, and B is what was parsed.
        parsed = os.stat(target)
        cached = module._DISK_QUARANTINE_CACHE[(str(tmp_path), self.LABEL)]
        assert cached[0][1] == (parsed.st_mtime_ns, parsed.st_ctime_ns, parsed.st_size,
                                parsed.st_ino), (
            "the answer was stored under the identity of a file that was never parsed: "
            "the signature came from one open and the content from another, so the key "
            "and the value beside it describe two different files")

    def test_restoring_the_original_is_a_cache_miss(self, monkeypatch, tmp_path):
        """
        The exploitation, measured rather than assumed -- and it does **not** reproduce at
        the base of this round.

        The brief expected restoring A to be a cache HIT returning B, and said `ctime`
        could not help because the cached one is A's. Measured, `ctime` is exactly what
        helps: `os.rename` moves the inode's `st_ctime_ns`, so the restored file no longer
        matches the stale key and the lookup misses. That is round 99's own repair paying
        for itself. The key describing a file that was never parsed is still wrong and is
        closed on its own terms above; this test pins the end-to-end answer so the claim
        does not drift in either direction.
        """
        family = tmp_path / "kinetics" / "families" / self.LABEL
        family.mkdir(parents=True)
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        a, b = self._bodies()
        manifest = family / QUARANTINE_FILENAME
        manifest.write_text(a)
        other = family / "manifest_b.py"
        other.write_text(b)

        state = {"fired": False}
        original = os.lstat
        target = os.path.abspath(str(manifest))

        def lstat(path, *args, **kwargs):
            answer = original(path, *args, **kwargs)
            if not state["fired"] and os.path.abspath(str(path)) == target:
                state["fired"] = True
                os.rename(target, str(family / "manifest_a.py"))
                os.rename(str(other), target)
            return answer

        monkeypatch.setattr(os, "lstat", lstat)
        resolve_quarantine(self.LABEL)
        monkeypatch.undo()
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        assert state["fired"] is True

        os.rename(target, str(other))
        os.rename(str(family / "manifest_a.py"), target)

        second = resolve_quarantine(self.LABEL)
        assert second[0] is not None and second[0].reason == "manifest A", (
            "restoring the original manifest returned the other file's quarantine")

    def test_the_cached_key_is_the_identity_of_the_file_that_was_parsed(
            self, monkeypatch, tmp_path):
        """
        The property directly, without the race: what is in the cache must describe the
        file whose bytes produced the value beside it.
        """
        from rmgpy.data.kinetics import quarantine as module

        family = tmp_path / "kinetics" / "families" / self.LABEL
        family.mkdir(parents=True)
        write_manifest(str(family))
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))

        assert resolve_quarantine(self.LABEL)[0] is not None
        cached = module._DISK_QUARANTINE_CACHE[(str(tmp_path), self.LABEL)]
        info = os.stat(str(family / QUARANTINE_FILENAME))
        assert cached[0][1] == (info.st_mtime_ns, info.st_ctime_ns, info.st_size,
                                info.st_ino)


class TestTheLoaderCarriesEveryFieldTheEntryHolds:
    """
    Round 102's HIGH 2. All three shapes in `get_library_reactions` named their fields by
    hand, so the template shape dropped `elementary_high_p`, `allow_pdep_route` and
    `allow_max_rate_violation`, and all three shapes dropped the last of those. A reaction
    that loses `elementary_high_p` silently misses pressure-dependent routing.

    The enumeration is derived from `Reaction` rather than written out, because a
    hand-kept list falling behind the class is the defect itself -- and a second
    hand-kept list in this file would reproduce it one layer up.
    """

    def _entry_declaring(self, **flags):
        """A library whose entry declares the flags, the way `load_entry` does."""
        return _library_declaring("a_seed", ["A_Family"], **flags)

    def test_the_template_shape_carries_the_three_flags(self):
        reaction = self._entry_declaring(
            elementary_high_p=True, allow_pdep_route=True,
            allow_max_rate_violation=True).get_library_reactions()[0]
        assert isinstance(reaction, TemplateReaction)
        assert reaction.elementary_high_p is True
        assert reaction.allow_pdep_route is True
        assert reaction.allow_max_rate_violation is True

    def test_the_ordinary_library_shape_carries_allow_max_rate_violation(self):
        """The field all three shapes dropped, on the shape that carried the other two."""
        reaction = self._entry_declaring(
            auto=False, allow_max_rate_violation=True).get_library_reactions()[0]
        assert isinstance(reaction, LibraryReaction)
        assert reaction.allow_max_rate_violation is True

    def test_a_false_flag_stays_false(self):
        """The negative control: carrying must not mean setting."""
        reaction = self._entry_declaring().get_library_reactions()[0]
        assert reaction.elementary_high_p is False
        assert reaction.allow_pdep_route is False
        assert reaction.allow_max_rate_violation is False

    def test_the_library_file_format_can_declare_all_three(self):
        """
        Through the real parser entry point, so the flags are shown to reach `entry.item`
        the way a shipped library file puts them there -- not only the way the fixture does.
        """
        library = KineticsLibrary(label="a_library", name="a_library")
        library.load_entry(index=1, label="r", kinetics=make_marcus(),
                           elementary_high_p=True, allow_pdep_route=True,
                           allow_max_rate_violation=True)
        item = library.entries[1].item
        assert item.elementary_high_p is True
        assert item.allow_pdep_route is True
        assert item.allow_max_rate_violation is True

    def test_every_field_of_reaction_is_carried_or_excluded_with_a_reason(self):
        """
        The partition is total, and machine-checked.

        A field added to `Reaction` tomorrow is carried by default. This asserts the
        classification covers everything, so a field can only be left behind by someone
        writing down why.
        """
        from rmgpy.data.kinetics import library as module

        carried = module.REACTION_STATE_FIELDS - set(module._NOT_CARRIED_FROM_ENTRY)
        excluded = set(module._NOT_CARRIED_FROM_ENTRY)

        assert len(module.REACTION_STATE_FIELDS) > 10, (
            "the field discovery found {0} fields on Reaction, which means the way it "
            "recognises them has stopped working and everything is now 'not a field' -- "
            "the silent-pass shape this partition exists to prevent".format(
                len(module.REACTION_STATE_FIELDS)))
        assert not excluded - module.REACTION_STATE_FIELDS, (
            "these names are excluded from the carry but are no longer fields of "
            "Reaction, so the exclusion list has gone stale: {0}".format(
                sorted(excluded - module.REACTION_STATE_FIELDS)))
        assert carried, "every field of Reaction is excluded; nothing is carried at all"

        unclassified = carried
        for name, reason in module._NOT_CARRIED_FROM_ENTRY.items():
            assert reason and len(reason) > 20, (
                "{0!r} is excluded without a reason worth reading".format(name))
        for expected in ("elementary_high_p", "allow_pdep_route",
                         "allow_max_rate_violation", "electrons", "degeneracy"):
            assert expected in unclassified, (
                "{0!r} must be carried from the entry".format(expected))


class TestCarryingTheDegeneracyDoesNotRestateTheRate:
    """
    Round 105's HIGH. `_carry_entry_fields` assigned `degeneracy` through the property,
    and `Reaction.degeneracy`'s setter is not an assignment: with kinetics already
    attached it multiplies the rate by a ratio and appends to the kinetics comment. The
    kinetics object at that moment **is** `entry.data`, the shared database object, so the
    edit reaches every later consumer of that entry and not only the reaction being built.

    Round 102's test used `degeneracy=3.0` on the template shape, where the ratio is
    exactly 1, and asserted the scalar only -- so it could see neither harm. The sweep is
    the repair to the test: the ratio has two branches, and which one fires depends on the
    value the constructor was given, which differs per shape.
    """

    A_BEFORE = 10.0
    SWEEP = (0.5, 1.0, 1.5, 2.0, 3.0)
    SHAPES = ("template", "originally_from", "ordinary")

    def _library(self, shape, degeneracy):
        """
        A one-entry library reaching each of the three branches of
        `get_library_reactions`, with the degeneracy declared on ``entry.item`` -- where
        `load_entry` puts it -- and an Arrhenius rate whose `A` is a round number.
        """
        if shape == "template":
            long_desc = ("Matched reaction 3 Lip + CH3 <=> CH3Li in A_Family/rate rule "
                         "[Root]\nEuclidian distance = 0\nfamily: A_Family")
        elif shape == "originally_from":
            long_desc = "Originally from reaction library: some_other_library"
        else:
            long_desc = ""
        library = KineticsLibrary(label="a_seed", name="a_seed")
        library.auto_generated = shape != "ordinary"
        library.entries = {
            1: Entry(
                index=1, label="Lip + CH3 <=> CH3Li",
                item=Reaction(
                    reactants=[Species(label="Lip", molecule=[Molecule(smiles="[Li+]")],
                                       reactive=False)],
                    products=[Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")],
                                      reactive=False)],
                    reversible=False, electrons=1, degeneracy=degeneracy),
                data=Arrhenius(A=(self.A_BEFORE, "m^3/(mol*s)"), n=0, Ea=(0, "kJ/mol"),
                               T0=(1, "K"), comment="Estimated from node Root"),
                long_desc=long_desc,
            )
        }
        return library

    @pytest.mark.parametrize("shape", SHAPES)
    @pytest.mark.parametrize("degeneracy", SWEEP)
    def test_the_carry_forwards_the_value_and_changes_nothing_else(self, shape,
                                                                   degeneracy):
        """
        The acceptance: the forwarded scalar, an unchanged `A`, and an unmutated entry.

        At `309acc0a4` the two `LibraryReaction` shapes pass no ``degeneracy`` to their
        constructor, so the old value is 1, the setter's ``< 2`` branch fires and the rate
        is multiplied by the **whole** new degeneracy -- 3.0 trebles it. On the template
        shape, which does pass it, only a degeneracy below 2 and not equal to 1 moves the
        rate. The comment is appended to in all fifteen combinations.
        """
        library = self._library(shape, degeneracy)
        entry = library.entries[1]
        comment_before = entry.data.comment

        reaction = library.get_library_reactions()[0]

        assert reaction.degeneracy == pytest.approx(degeneracy), (
            "the declared degeneracy did not reach the reaction at all")
        assert reaction.kinetics.A.value_si == pytest.approx(self.A_BEFORE), (
            "carrying degeneracy={0} onto the {1} shape rescaled the rate by {2}; a "
            "library file that declares a degeneracy is stating the rate it wants, and "
            "re-deriving the rate from a number that was already true applies it "
            "twice".format(degeneracy, shape,
                           reaction.kinetics.A.value_si / self.A_BEFORE))
        assert entry.data.A.value_si == pytest.approx(self.A_BEFORE), (
            "the rate was changed in place on entry.data -- the shared database object -- "
            "so every later consumer of this entry sees the corrupted value, not just "
            "the reaction being loaded")
        assert entry.data.comment == comment_before, (
            "the setter appended to the shared kinetics comment. This campaign reads "
            "authorship out of that comment, so it is not cosmetic")

    def test_the_rate_is_shared_with_the_entry_so_the_damage_would_escape(self):
        """
        The premise the two assertions above rest on, measured rather than assumed.

        If the loader handed the reaction a *copy* of `entry.data`, mutating it would be
        confined to the reaction and the second assertion would be vacuous.
        """
        library = self._library("ordinary", 2.0)
        entry = library.entries[1]
        reaction = library.get_library_reactions()[0]
        assert reaction.kinetics is entry.data, (
            "the reaction no longer shares the entry's kinetics object, so the test above "
            "measures nothing; either the loader started copying or the fixture broke")

    def test_the_storage_table_names_a_setter_that_really_transforms(self):
        """
        The negative control on the workaround, so it cannot outlive its reason.

        `_CARRIED_THROUGH_STORAGE` exists only because assigning `degeneracy` through the
        property transforms the rate. If upstream ever makes that setter a plain
        assignment, this fails and says to delete the special case rather than carrying an
        unexplained back door forever.

        **Structural** at `309acc0a4`, which has no `_CARRIED_THROUGH_STORAGE`. A guard on
        the repair rather than evidence of the defect; the fifteen sweep cases in this
        class carry the behavioural red.
        """
        from rmgpy.data.kinetics import library as module

        for name, storage in module._CARRIED_THROUGH_STORAGE.items():
            assert name in module.REACTION_STATE_FIELDS, (
                "{0!r} is routed through storage but is no longer a field of "
                "Reaction".format(name))
            assert hasattr(Reaction(), storage), (
                "{0!r} is routed to {1!r}, which Reaction no longer has".format(
                    name, storage))

        probe = Reaction(kinetics=Arrhenius(A=(self.A_BEFORE, "m^3/(mol*s)"), n=0,
                                            Ea=(0, "kJ/mol"), T0=(1, "K")),
                         degeneracy=1)
        probe.degeneracy = 4.0
        assert probe.kinetics.A.value_si != pytest.approx(self.A_BEFORE), (
            "assigning degeneracy through the property no longer changes the rate, so "
            "_CARRIED_THROUGH_STORAGE's degeneracy entry has no reason left to exist")


class TestTheConversionEnumeratesTheSourceNotTheDestination:
    """
    Round 105's first MEDIUM. `as_library_reaction` enumerated
    `LibraryReaction.__init__`'s parameters, so state the class holds and the constructor
    does not take was dropped by construction -- `is_forward`, `rank`, `comment` and
    `label`. Round 102's test iterated the same constructor signature, so the enumeration
    and its check shared a blind spot exactly.

    The enumeration comes from the source's class now, and this file checks it against
    that class rather than against a second list of names.
    """

    def _markers(self):
        """
        One distinctive value per carried field, checked against the class below.

        A name here that `Reaction` no longer has, or a carried field with no marker,
        fails `test_every_carried_field_arrives`. That is what makes a field added to
        `Reaction` tomorrow carried **or loudly refused**, rather than quietly dropped.
        """
        collider = Species(label="Ar", molecule=[Molecule(smiles="[Ar]")])
        return {
            "allow_max_rate_violation": True,
            "allow_pdep_route": True,
            "comment": "a comment worth keeping",
            "degeneracy": 2.0,
            "duplicate": True,
            "electrons": -1,
            "elementary_high_p": True,
            "index": 11,
            "is_forward": True,
            "kinetics": Arrhenius(A=(10.0, "m^3/(mol*s)"), n=0, Ea=(0, "kJ/mol"),
                                  T0=(1, "K"), comment="the rate"),
            "label": "Lip + CH3 <=> CH3Li",
            "network_kinetics": Arrhenius(A=(3.0, "m^3/(mol*s)"), n=0, Ea=(0, "kJ/mol"),
                                          T0=(1, "K")),
            "pairs": [("Lip", "CH3Li")],
            "rank": 7,
            "reversible": False,
            "specific_collider": collider,
            "transition_state": TransitionState(),
        }

    def _source(self, markers):
        """
        A template reaction holding every marker, built the way production does.

        `kinetics` is attached **last** on purpose: assigning `degeneracy` while kinetics
        are attached is the very transformation this round is about, and a fixture that
        triggered it would be measuring itself.
        """
        source = TemplateReaction(
            reactants=[Species(label="Lip", molecule=[Molecule(smiles="[Li+]")])],
            products=[Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")])],
            family="A_Family", template="Root")
        for name, value in markers.items():
            if name == "kinetics":
                continue
            setattr(source, name, value)
        source.kinetics = markers["kinetics"]
        return source

    def test_the_four_fields_the_constructor_does_not_take_survive(self):
        """The named acceptance, on its own, so a regression says which fields went."""
        from rmgpy.rmg.model import as_library_reaction

        markers = self._markers()
        converted = as_library_reaction(self._source(markers), "a_library")
        for name in ("is_forward", "rank", "comment", "label"):
            assert getattr(converted, name) == markers[name], (
                "{0!r} is state LibraryReaction holds and its constructor does not take, "
                "so enumerating the constructor dropped it: {1!r} -> {2!r}".format(
                    name, markers[name], getattr(converted, name)))

    def test_every_carried_field_arrives(self):
        """
        The general form: every field the partition calls carried, checked by value.

        Iterating `REACTION_STATE_FIELDS` rather than a list of names is the point -- a
        field added to `Reaction` tomorrow enters this test automatically and fails it
        until it is given a marker and shown to arrive.

        **Structural** at `309acc0a4`, which has no `REACTION_STATE_FIELDS`: this pins
        the shape of the enumeration, and
        `test_the_four_fields_the_constructor_does_not_take_survive` carries the
        behavioural red for the same defect.
        """
        from rmgpy.data.kinetics import library as library_module
        from rmgpy.rmg import model as model_module
        from rmgpy.rmg.model import as_library_reaction

        carried = (library_module.REACTION_STATE_FIELDS
                   - set(model_module._NOT_CARRIED_IN_CONVERSION))
        markers = self._markers()
        assert set(markers) == carried, (
            "the markers and the carried set have diverged. missing markers: {0}; "
            "markers for fields that are no longer carried: {1}".format(
                sorted(carried - set(markers)), sorted(set(markers) - carried)))

        converted = as_library_reaction(self._source(markers), "a_library")
        for name in sorted(carried):
            assert getattr(converted, name) == markers[name], (
                "{0!r} did not survive the conversion: {1!r} -> {2!r}".format(
                    name, markers[name], getattr(converted, name)))

    def test_the_conversion_partition_is_total_and_reasoned(self):
        """
        **Structural** at `309acc0a4`: it names `_NOT_CARRIED_IN_CONVERSION`, which the
        base does not have, so it fails there with an `AttributeError` rather than on
        the defect. It pins the partition's shape; the behavioural evidence for the
        same ground is `test_the_four_fields_the_constructor_does_not_take_survive`.
        """
        from rmgpy.data.kinetics import library as library_module
        from rmgpy.rmg import model as model_module

        excluded = set(model_module._NOT_CARRIED_IN_CONVERSION)
        assert not excluded - library_module.REACTION_STATE_FIELDS, (
            "these names are excluded from the conversion but are no longer fields of "
            "Reaction: {0}".format(
                sorted(excluded - library_module.REACTION_STATE_FIELDS)))
        for name, reason in model_module._NOT_CARRIED_IN_CONVERSION.items():
            assert reason and len(reason) > 20, (
                "{0!r} is excluded from the conversion without a reason worth "
                "reading".format(name))

    # `test_every_field_called_carried_can_actually_be_assigned` stood here until round
    # 107. It probed a bare `Reaction` against itself, which is not a pair any carry runs
    # on, so a descriptor only the subclasses disagreed about would have passed it. Its
    # replacement iterates the real (source, target) class pairs and runs the real carry:
    # `TestACarryThatCannotCarryRefusesLoudly`.

    def test_the_conversion_does_not_rescale_the_rate(self):
        """
        The guard on the repair itself: the conversion now carries `degeneracy` too, and
        it must carry it the same way the loader does. Green at `309acc0a4` -- the
        constructor path never rescaled -- and red the moment the storage rule is dropped.
        """
        from rmgpy.rmg.model import as_library_reaction

        markers = self._markers()
        markers["degeneracy"] = 0.5
        source = self._source(markers)
        before = source.kinetics.A.value_si
        converted = as_library_reaction(source, "a_library")
        assert converted.kinetics.A.value_si == pytest.approx(before)
        assert converted.degeneracy == pytest.approx(0.5)

    def test_the_conversion_still_replaces_the_library_slot(self):
        """The negative control: carrying everything must not carry the family label."""
        from rmgpy.rmg.model import as_library_reaction

        converted = as_library_reaction(self._source(self._markers()), "a_library")
        assert isinstance(converted, LibraryReaction)
        assert converted.library == "a_library"
        assert converted.family == "a_library", (
            "`family` means the LIBRARY on this class, which is what the conversion is "
            "for; carrying the source's family label here would undo it")

    def test_the_conversion_does_not_share_the_reactant_list(self):
        from rmgpy.rmg.model import as_library_reaction

        source = self._source(self._markers())
        converted = as_library_reaction(source, "a_library")
        assert converted.reactants is not source.reactants
        assert converted.products is not source.products


class TestTheAbsenceOfONofollowIsRefused:
    """
    Round 105's second MEDIUM. `getattr(os, 'O_NOFOLLOW', 0)` degrades to a flag of 0,
    which is not a weaker containment but none: every component below the anchor would be
    followed and a manifest outside the database read and executed, with nothing said.
    Round 101 enforced the `dir_fd` half of the same platform requirement and left this
    half documented -- and a documented hole is still a hole.
    """

    def setup_method(self):
        _clear_gate_caches()

    def teardown_method(self):
        _clear_gate_caches()

    def test_a_platform_without_o_nofollow_is_refused(self, monkeypatch, tmp_path):
        label = self._escape_at(monkeypatch, tmp_path, "families")
        monkeypatch.delattr(os, "O_NOFOLLOW")

        assert resolve_quarantine(label) == (None, False), (
            "with O_NOFOLLOW absent the flag became 0, the descent followed the "
            "symlinked `kinetics/families` and the manifest behind it was executed; "
            "there is no containment available on that path, so it must refuse rather "
            "than answer")

    def test_the_same_escape_is_refused_with_the_flag_present(self, monkeypatch,
                                                              tmp_path):
        """The control: the refusal above must not be the only thing being measured."""
        label = self._escape_at(monkeypatch, tmp_path, "families")
        assert resolve_quarantine(label) == (None, False)

    def test_an_in_tree_family_is_still_answered(self, monkeypatch, tmp_path):
        """And the flag being present must still let an ordinary family through."""
        _quarantined_on_disk(tmp_path, QUARANTINED_FIRST)
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        quarantine, answered = resolve_quarantine(QUARANTINED_FIRST)
        assert answered and quarantine is not None

    _escape_at = TestEveryComponentBelowTheAnchorIsPinned._escape_at


#: How to compare a field before and after a transform. The default is ``==``; a field
#: appears here only because the object it holds has no value equality, so a deep copy of
#: it would never compare equal to its original however faithfully it was reproduced.
#: A field added to `Reaction` tomorrow gets ``==`` without anyone deciding to give it
#: one, and fails loudly if that is the wrong answer -- which is the point.
_COMPARED_BY = {
    "reactants": lambda v: [s.label for s in v],
    "products": lambda v: [s.label for s in v],
    "specific_collider": lambda v: None if v is None else v.label,
    "kinetics": lambda v: None if v is None else (v.A.value_si, v.n.value_si),
    "network_kinetics": lambda v: None if v is None else (v.A.value_si, v.n.value_si),
    "transition_state": lambda v: type(v).__name__,
    "entry": lambda v: None if v is None else (v.index, v.label),
    "reverse": lambda v: None if v is None else v.index,
    "labeled_atoms": lambda v: sorted(v),
    # By value, so a *detached* copy still compares equal here. Identity is asserted
    # separately, by `TestTheCopyIsOneGraphNotFourLists`, because it is exactly the
    # property a value comparison structurally cannot see.
    "pairs": lambda v: [tuple(s.label for s in pair) for pair in v],
}


def _compare(name, value):
    return _COMPARED_BY.get(name, lambda v: v)(value)


class TestEveryTransformReproducesTheWholeState:
    """
    Round 107's HIGH. Rounds 102 and 105 made the *loader* preserve every field a
    `Reaction` carries by discovering the field set from the class. Two transforms still
    enumerated by hand, so the preservation ended at the first ``pickle`` or ``copy()``:
    measured at `5a821898d`, `TemplateReaction.__reduce__` dropped `elementary_high_p`,
    `allow_pdep_route`, `allow_max_rate_violation`, `rank`, `comment` and `label`, and
    `copy()` dropped the same three flags plus `rank`, `network_kinetics` and
    `labeled_atoms` while leaving `comment` as ``None`` where the class declares a
    ``str``.

    The census is wider than those two methods. `LibraryReaction.__reduce__` dropped
    `rank`, `comment`, `label` and turned `is_forward` ``True`` into ``False`` -- round
    105's three fields, one transform over -- and `LibraryReaction` had no `copy()` at
    all, so it inherited `Reaction.copy` and returned a **base `Reaction`**, losing
    `library`, `family` and `entry`, the carrier the gate reads authorship from.

    All four go through one derived definition of what reaction state is now, so the next
    field cannot be added to one site and forgotten in three. These tests iterate that
    definition rather than a list of names, so a field added to `Reaction` tomorrow
    enters them automatically.
    """

    SHAPES = ("template", "library")
    TRANSFORMS = ("pickle", "copy")

    def _markers(self):
        """One distinctive value per reproduced field, checked against the class below."""
        collider = Species(label="Ar", molecule=[Molecule(smiles="[Ar]")])
        return {
            "allow_max_rate_violation": True,
            "allow_pdep_route": True,
            "comment": "a comment worth keeping",
            "degeneracy": 3.0,
            "duplicate": True,
            "electrons": -1,
            "elementary_high_p": True,
            "index": 11,
            "is_forward": True,
            "kinetics": Arrhenius(A=(10.0, "m^3/(mol*s)"), n=0, Ea=(0, "kJ/mol"),
                                  T0=(1, "K"), comment="the rate"),
            "label": "Lip + CH3 <=> CH3Li",
            "network_kinetics": Arrhenius(A=(3.0, "m^3/(mol*s)"), n=0, Ea=(0, "kJ/mol"),
                                          T0=(1, "K")),
            # `pairs` and `labeled_atoms` are filled in by `_built`, from the species it
            # owns. String stand-ins stood here until round 108 and were exactly why the
            # tests could not see a severed cross-reference: a string cannot point into a
            # copied structure, so no assertion over one can notice that it stopped.
            "pairs": None,
            "products": None,          # filled in by `_built`, which owns the species
            "rank": 7,
            "reactants": None,         # likewise
            "reversible": False,
            "specific_collider": collider,
            "transition_state": TransitionState(),
        }

    def _built(self, shape):
        """
        A fully populated reaction of either shape, plus the markers it was given.

        `kinetics` is attached **last**: assigning `degeneracy` while kinetics are
        attached is the transformation round 105 was about, and a fixture that triggered
        it would be measuring itself.
        """
        markers = self._markers()
        reactants = [Species(label="Lip", molecule=[Molecule(smiles="[Li+]")])]
        products = [Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")])]
        markers["reactants"] = reactants
        markers["products"] = products
        markers["pairs"] = [(reactants[0], products[0])]
        labelled = reactants[0].molecule[0].atoms[0]
        labelled.label = "*1"
        entry = Entry(index=4, label="an entry", long_desc="family: A_Family")

        if shape == "template":
            reaction = TemplateReaction(reactants=reactants, products=products,
                                        family="A_Family", entry=entry)
            markers.update({
                "family": "A_Family",
                "template": ["Root"],
                "estimator": "rate rules",
                "reverse": TemplateReaction(index=99),
                "entry": entry,
                "labeled_atoms": {"reactants": {"*1": labelled}, "products": {}},
            })
        else:
            reaction = LibraryReaction(reactants=reactants, products=products,
                                       library="a_library", entry=entry)
            markers.update({
                "family": "a_library",
                "library": "a_library",
                "entry": entry,
            })

        for name, value in markers.items():
            if name == "kinetics":
                continue
            setattr(reaction, name, value)
        reaction.kinetics = markers["kinetics"]
        return reaction, markers

    def _transformed(self, transform, reaction):
        if transform == "pickle":
            return pickle.loads(pickle.dumps(reaction))
        return reaction.copy()

    def _reproduced(self, reaction):
        from rmgpy.data.kinetics.family import _NOT_REPRODUCED, state_fields

        return state_fields(reaction) - set(_NOT_REPRODUCED)

    @pytest.mark.parametrize("transform", TRANSFORMS)
    @pytest.mark.parametrize("shape", SHAPES)
    def test_the_three_flags_and_the_carrier_survive(self, shape, transform):
        """
        The named acceptance, on its own, so a regression says which field went.

        `elementary_high_p` is the one with teeth: a reaction that loses it misses
        pressure-dependent routing, and nothing says so. `is_forward`, `rank` and
        `comment` are here because they are round 105's three fields, which the library
        shape's pickle was dropping all over again.

        This test imports nothing the repair adds, so it runs at the base too and its
        red state is the dropped field rather than a missing constant.
        """
        reaction, markers = self._built(shape)
        after = self._transformed(transform, reaction)

        for name in ("elementary_high_p", "allow_pdep_route", "allow_max_rate_violation",
                     "degeneracy", "is_forward", "rank", "comment", "label"):
            assert getattr(after, name) == markers[name], (
                "{0!r} did not survive {1} of a {2} reaction: {3!r} -> {4!r}".format(
                    name, transform, shape, markers[name], getattr(after, name)))
        assert after.entry is not None and after.entry.label == "an entry", (
            "the entry carrier did not survive {0} of a {1} reaction".format(
                transform, shape))
        assert type(after) is type(reaction), (
            "{0} of a {1} reaction returned a {2}".format(
                transform, shape, type(after).__name__))

    @pytest.mark.parametrize("transform", TRANSFORMS)
    @pytest.mark.parametrize("shape", SHAPES)
    def test_every_reproduced_field_arrives(self, shape, transform):
        """
        The general form, and the acceptance for "a newly-added field cannot be dropped".

        The field set is discovered from the object, not listed here, so a field added to
        `Reaction` tomorrow is checked by this test without anyone remembering to add it
        -- and fails it until it is given a marker and shown to survive both transforms.

        **Structural** at `5a821898d`, which has no `state_fields` or `_NOT_REPRODUCED`:
        it pins the shape of the definition. `test_the_three_flags_and_the_carrier_survive`
        carries the behavioural red for the same four combinations, naming every field
        that was dropped.
        """
        reaction, markers = self._built(shape)
        reproduced = self._reproduced(reaction)

        assert reproduced == set(markers), (
            "the markers and the reproduced state have diverged. fields with no marker: "
            "{0}; markers for fields no longer reproduced: {1}".format(
                sorted(reproduced - set(markers)), sorted(set(markers) - reproduced)))

        after = self._transformed(transform, reaction)
        for name in sorted(reproduced):
            assert _compare(name, getattr(after, name)) == _compare(name, markers[name]), (
                "{0!r} did not survive {1} of a {2} reaction: {3!r} -> {4!r}".format(
                    name, transform, shape, markers[name], getattr(after, name)))

    @pytest.mark.parametrize("shape", SHAPES)
    def test_the_copy_is_deep_where_it_says_it_is(self, shape):
        """
        The other half of `copy()`: carried-by-reference is right for a label and wrong
        for a list. Round 105's HIGH escaped into the database through a shared kinetics
        object, so a copy that shared one would be the same defect wearing a new hat.

        **Structural** at `5a821898d`, which has neither exclusion table. The base's
        `copy()` did deepen the fields it copied at all; what it did instead was drop
        `network_kinetics` and `labeled_atoms` entirely, which the probe measures.
        """
        from rmgpy.data.kinetics.family import (_NOT_REPRODUCED,
                                                _TEMPLATE_NOT_COPIED_BY_REFERENCE)

        reaction, _ = self._built(shape)
        other = reaction.copy()
        for name in sorted(set(_TEMPLATE_NOT_COPIED_BY_REFERENCE) - set(_NOT_REPRODUCED)):
            mine = getattr(reaction, name, None)
            if mine is None:
                continue
            assert mine is not getattr(other, name), (
                "copy() shares {0!r} with the original, so an edit to one reaches the "
                "other".format(name))

    @pytest.mark.parametrize("shape", SHAPES)
    def test_the_copy_does_not_restate_the_rate(self, shape):
        """
        The guard on the repair itself. `copy()` carries `degeneracy` through the generic
        helper now, and the helper must write the storage: assigning the property while
        kinetics are attached would multiply the copied rate by the degeneracy all over
        again. Green before the repair because the old `copy()` assigned degeneracy while
        `kinetics` was still unset -- safety that was positional and unmarked.
        """
        reaction, markers = self._built(shape)
        other = reaction.copy()
        assert other.degeneracy == 3.0
        assert other.kinetics.A.value_si == pytest.approx(
            markers["kinetics"].A.value_si), (
            "copy() restated the rate while carrying the degeneracy")

    def test_the_sites_that_enumerate_reaction_state_share_one_definition(self):
        """
        The census, as a check rather than a claim.

        Five sites in the gated files enumerate reaction state: the loader's
        `_carry_entry_fields`, and `__reduce__` and `copy` on each of the two reaction
        classes. Each must reach the shared helpers, and none may read ``self.<field>``
        for a field the helper is responsible for -- which is what a re-added hand
        enumeration looks like.

        This is syntactic. It cannot tell that a site calls the helper with the *right*
        policy, only that it does not enumerate; the value checks above cover the rest.

        **Structural** at `5a821898d` by construction -- it is the census of the repair.
        The behavioural consequence of the census being wrong is what every other test
        in this class measures.
        """
        from rmgpy.data.kinetics import family as family_module
        from rmgpy.data.kinetics import library as library_module

        # `copy_reaction` counts because it *is* the shared definition for the copies --
        # it reaches the other three itself. Round 108 moved both `copy()` bodies into it
        # and this check went red on the delegation, which was the checker being too
        # literal rather than a site enumerating by hand again.
        shared = ("reaction_state", "carry_reaction_state", "apply_reaction_state",
                  "copy_reaction")
        sites = {
            "TemplateReaction.__reduce__": TemplateReaction.__reduce__,
            "TemplateReaction.copy": TemplateReaction.copy,
            "LibraryReaction.__reduce__": LibraryReaction.__reduce__,
            "LibraryReaction.copy": LibraryReaction.copy,
            "_carry_entry_fields": library_module._carry_entry_fields,
        }
        assert len(sites) == 5, "the census has changed size without the count moving"

        deepened = set(family_module._TEMPLATE_NOT_COPIED_BY_REFERENCE) - set(
            family_module._NOT_REPRODUCED)
        for name, site in sites.items():
            source = inspect.getsource(site)
            body = source.split('"""')[-1]
            assert any(helper in source for helper in shared), (
                "{0} does not reach the shared definition of reaction state, so it is "
                "enumerating by hand again".format(name))
            # On a word boundary: a plain substring test reads `self.labeled_atoms` as a
            # hand-enumeration of `label`, which is the checker crying wolf rather than
            # the defect it was written for.
            enumerated = sorted(
                field for field in family_module.REACTION_STATE_FIELDS
                if re.search(r"self\.{0}\b".format(re.escape(field)), body)
                and field not in deepened)
            assert not enumerated, (
                "{0} reads {1} off the object by hand; that is the enumeration this "
                "round removed".format(name, enumerated))

    def test_the_definition_is_discovered_from_both_halves_of_the_class(self):
        """
        `state_fields` must not quietly stop finding things. A discovery that returned an
        empty set would make every test above pass vacuously.

        **Structural** at `5a821898d`: `state_fields` is this round's function. It is the
        anti-vacuity guard on the tests above, not evidence of the defect.
        """
        from rmgpy.data.kinetics.family import REACTION_STATE_FIELDS, state_fields

        assert len(REACTION_STATE_FIELDS) >= 20, (
            "the class-half of the discovery has stopped finding fields")
        reaction, _ = self._built("template")
        found = state_fields(reaction)
        assert {"family", "template", "estimator", "reverse", "entry",
                "labeled_atoms"} <= found, (
            "the instance-half of the discovery misses TemplateReaction's own fields")
        assert state_fields(Reaction()) == REACTION_STATE_FIELDS, (
            "a base Reaction has no instance dictionary, so the two halves must agree")


class TestACarryThatCannotCarryRefusesLoudly:
    """
    Round 107's MEDIUM. `carry_reaction_state` caught `AttributeError` from both the
    source read and the target write and continued, so a field the partition called
    carried could be dropped in silence -- the opposite of what its docstring promised,
    and exactly the shape (`protons`, round 102) the docstring was written about.

    It raises now. The partition test below is what makes that safe: no field any policy
    calls carried can reach the raise, so the loud failure costs nothing a correct
    partition does not already forbid, and announces the partition going stale the moment
    it stops being correct.
    """

    #: The (source, target) class pairs the carries actually run on in production. Round
    #: 105's version of this test probed a bare `Reaction` against itself, so a
    #: descriptor only the subclasses disagreed about would have passed it.
    PAIRS = (
        (Reaction, TemplateReaction),          # the loader, template shape
        (Reaction, LibraryReaction),           # the loader, both library shapes
        (TemplateReaction, LibraryReaction),   # rmg/model.py's conversion
        (TemplateReaction, TemplateReaction),  # __reduce__ / copy
        (LibraryReaction, LibraryReaction),    # __reduce__ / copy
    )

    def _policies(self):
        from rmgpy.data.kinetics import family as family_module
        from rmgpy.data.kinetics import library as library_module
        from rmgpy.rmg import model as model_module

        return {
            "the loader": library_module._NOT_CARRIED_FROM_ENTRY,
            "the conversion": model_module._NOT_CARRIED_IN_CONVERSION,
            "the transforms": family_module._NOT_REPRODUCED,
            "copy()": family_module._TEMPLATE_NOT_COPIED_BY_REFERENCE,
        }

    def test_every_field_called_carried_can_be_carried_between_the_real_classes(self):
        """
        Carried must mean carried, on the classes the carry really runs between.

        **Structural** at `5a821898d`, which has no `_NOT_REPRODUCED`. It would also
        have *passed* there if it could import: the base swallowed the refusal instead
        of raising, so silence looked like success. That is the MEDIUM in one sentence,
        and `test_a_field_the_target_refuses_is_refused_loudly` is what measures it.
        """
        from rmgpy.data.kinetics.family import REACTION_STATE_FIELDS, carry_reaction_state

        for policy in self._policies().values():
            for source_class, target_class in self.PAIRS:
                source, target = source_class(), target_class()
                for name in sorted(REACTION_STATE_FIELDS - set(policy)):
                    carry_reaction_state(target, source, policy, fields={name})

    def test_a_field_the_target_refuses_is_refused_loudly(self):
        """
        The acceptance: take the read-only field out of a real policy and the carry must
        say so rather than continue. `protons` is read-only because it is derived from
        the charge balance of the reactants and products.

        Written to run **at the base as well as here** -- it imports only names that
        already existed and names the exception by string -- so its red state is the
        silence itself (``DID NOT RAISE``) rather than an `ImportError` for a constant
        the repair adds.
        """
        from rmgpy.data.kinetics.library import (_NOT_CARRIED_FROM_ENTRY,
                                                 carry_reaction_state)

        claimed = {name: reason for name, reason in _NOT_CARRIED_FROM_ENTRY.items()
                   if name != "protons"}
        with pytest.raises(Exception) as raised:
            carry_reaction_state(LibraryReaction(), Reaction(), claimed)
        assert type(raised.value).__name__ == "ReactionStateNotCarried"
        assert "protons" in str(raised.value)

    def test_a_source_that_holds_none_of_the_state_is_refused_loudly(self):
        """
        The other direction, and the one that was truly unreachable before: with the
        source read swallowed too, a carry from an object holding none of the state
        completed in silence and produced a reaction with nothing on it.
        """
        from rmgpy.data.kinetics.library import (_NOT_CARRIED_FROM_ENTRY,
                                                 carry_reaction_state)

        class HoldsNoReactionState:
            pass

        with pytest.raises(Exception) as raised:
            carry_reaction_state(LibraryReaction(), HoldsNoReactionState(),
                                 _NOT_CARRIED_FROM_ENTRY)
        assert type(raised.value).__name__ == "ReactionStateNotCarried"

    def test_the_refusal_cannot_be_swallowed_by_an_attributeerror_handler(self):
        """
        The refusal must not be an `AttributeError`. Both failures it reports *are*
        `AttributeError`s underneath, and the call sites sit inside loaders full of
        ``except AttributeError`` -- an exception one of those could catch would restore
        the silence this round removed.
        """
        from rmgpy.data.kinetics.library import (_NOT_CARRIED_FROM_ENTRY,
                                                 carry_reaction_state)

        claimed = {name: reason for name, reason in _NOT_CARRIED_FROM_ENTRY.items()
                   if name != "protons"}
        try:
            carry_reaction_state(LibraryReaction(), Reaction(), claimed)
        except AttributeError as error:
            pytest.fail("the refusal is an AttributeError, which the loaders around it "
                        "catch and discard: {0!r}".format(error))
        except Exception:
            return
        pytest.fail("a field the target refuses was carried in silence")

    def test_the_documentation_says_it_refuses(self):
        """
        The code and the docstring must agree; they did not before this round, which is
        the MEDIUM in one line. Imported from `library`, where the function was when the
        docstring made the opposite promise.
        """
        from rmgpy.data.kinetics import library as library_module

        checked = 0
        for name in ("carry_reaction_state", "reaction_state", "apply_reaction_state"):
            function = getattr(library_module, name, None)
            if function is None:
                continue
            checked += 1
            assert "raises" in (function.__doc__ or "").lower(), (
                "{0} refuses a field it cannot carry and does not say so".format(name))
        assert checked, "none of the carrying helpers is reachable under its own name"


class TestTheCopyIsOneGraphNotFourLists:
    """
    Round 108's HIGH. A reaction is not a bag of independent fields. `pairs` holds the
    reaction's own `Species` objects -- `Reaction.generate_pairs` appends them straight out
    of `self.reactants` and `self.products` -- and `labeled_atoms` holds `Atom` objects out
    of those species' molecules. Its internal consistency is a property of the references
    *between* those four lists, not of the lists separately.

    `copy()` deep-copied each of them in its own call, each with its own memo, which severs
    exactly those references. Measured at `56156ac9a`: after a copy, `pairs[0][0]` was in
    neither the copy's reactants nor the original's, and relabelling through
    `labeled_atoms` changed nothing any structure could see. `Species.__eq__` is identity,
    so `reactants.index(pair[0])` raises on the first; the second is silent.

    Every assertion here is `is`, deliberately. Round 107's tests compare by value and are
    green at `56156ac9a` throughout -- a detached clone has all the right values. Value
    equality is structurally unable to see this defect, which is the second half of the
    finding and the reason these fixtures hold real species and real atoms rather than the
    strings they held until this round.
    """

    TRANSFORMS = ("copy", "pickle")

    def _built(self, shape):
        """A reaction whose `pairs` and `labeled_atoms` point into its own species."""
        reactants = [Species(label="ethane", molecule=[Molecule(smiles="CC")])]
        products = [Species(label="ethyl", molecule=[Molecule(smiles="C[CH2]")])]
        if shape == "template":
            reaction = TemplateReaction(reactants=reactants, products=products,
                                        family="A_Family")
        else:
            reaction = LibraryReaction(reactants=reactants, products=products,
                                       library="a_library")
        reaction.pairs = [(reactants[0], products[0])]
        labelled = reactants[0].molecule[0].atoms[0]
        labelled.label = "*1"
        if shape == "template":
            reaction.labeled_atoms = {"reactants": {"*1": labelled}, "products": {}}
        return reaction

    def _transformed(self, transform, reaction):
        if transform == "pickle":
            return pickle.loads(pickle.dumps(reaction))
        return reaction.copy()

    @staticmethod
    def _atoms_of(reaction):
        return [atom
                for species in list(reaction.reactants) + list(reaction.products)
                for molecule in species.molecule
                for atom in molecule.atoms]

    def test_the_fixture_holds_the_invariant_before_anything_is_copied(self):
        """
        Anti-vacuity. If the fixture did not point into its own species -- as the string
        version did not -- every assertion below would be about nothing.

        **Green at `56156ac9a`** by construction: it is the guard on the other tests, not
        evidence of the defect.
        """
        for shape in ("template", "library"):
            reaction = self._built(shape)
            assert any(reaction.pairs[0][0] is s for s in reaction.reactants)
            assert any(reaction.pairs[0][1] is s for s in reaction.products)
            if shape == "template":
                assert any(reaction.labeled_atoms["reactants"]["*1"] is a
                           for a in self._atoms_of(reaction))

    @pytest.mark.parametrize("transform", TRANSFORMS)
    @pytest.mark.parametrize("shape", ("template", "library"))
    def test_every_pair_member_is_one_of_the_copys_own_species(self, shape, transform):
        """Every `pairs` member must be a species the transformed reaction owns.

        The ``copy`` arms are the behavioural red state at `56156ac9a`. The ``pickle``
        arms are **green at the base**: one state dict, one pickle memo, so the
        references survive by construction and always did. Pinned, not claimed.
        """
        reaction = self._built(shape)
        after = self._transformed(transform, reaction)
        owned = list(after.reactants) + list(after.products)
        original = list(reaction.reactants) + list(reaction.products)
        for pair in after.pairs:
            for member in pair:
                assert any(member is s for s in owned), (
                    "after {0} of a {1} reaction, a pairs member is a species the "
                    "reaction does not own; Species.__eq__ is identity, so "
                    "reactants.index(pair[0]) raises".format(transform, shape))
                assert not any(member is s for s in original), (
                    "after {0}, a pairs member is still the *original* reaction's "
                    "species, so the copy is not a copy".format(transform))

    @pytest.mark.parametrize("transform", TRANSFORMS)
    def test_every_labelled_atom_is_inside_one_of_the_copys_own_species(self, transform):
        """Every labelled atom must be an atom of one of the copy's own molecules.

        The ``copy`` arms are the behavioural red state at `56156ac9a`. The ``pickle``
        arms are **green at the base**: one state dict, one pickle memo, so the
        references survive by construction and always did. Pinned, not claimed.
        """
        reaction = self._built("template")
        after = self._transformed(transform, reaction)
        owned = self._atoms_of(after)
        original = self._atoms_of(reaction)
        for group in after.labeled_atoms.values():
            for label, atom in group.items():
                for one in (atom if isinstance(atom, list) else [atom]):
                    assert any(one is a for a in owned), (
                        "after {0}, labelled atom {1!r} is not inside any species the "
                        "copy owns, so relabelling through it reaches nothing".format(
                            transform, label))
                    assert not any(one is a for a in original), (
                        "after {0}, labelled atom {1!r} is still the original's "
                        "atom".format(transform, label))

    @pytest.mark.parametrize("transform", TRANSFORMS)
    def test_relabelling_through_labeled_atoms_is_observable_in_the_structure(
            self, transform):
        """
        The consequence rather than the mechanism. `family.py:2593` relabels the reaction's
        structures out of `labeled_atoms` before regenerating `pairs` and `template`; a
        detached atom makes that a no-op that reports nothing.

        The ``copy`` arms are the behavioural red state at `56156ac9a`. The ``pickle``
        arms are **green at the base**: one state dict, one pickle memo, so the
        references survive by construction and always did. Pinned, not claimed.
        """
        reaction = self._built("template")
        after = self._transformed(transform, reaction)
        after.labeled_atoms["reactants"]["*1"].label = "*9"

        assert any(a.label == "*9" for a in self._atoms_of(after)), (
            "relabelling through labeled_atoms after a {0} changed nothing the reaction "
            "can see".format(transform))
        assert all(a.label != "*9" for a in self._atoms_of(reaction)), (
            "relabelling the copy reached back into the original")

    @pytest.mark.parametrize("shape", ("template", "library"))
    def test_the_shallow_half_is_still_shallow(self, shape):
        """
        The repair must not deepen what was deliberately shallow. `entry` is the shared
        database object the quarantine gate reads authorship from; a copy that cloned it
        would be a different defect in the same place.

        **Green at `56156ac9a`** -- a guard against the repair overshooting.
        """
        entry = Entry(index=4, label="an entry", long_desc="family: A_Family")
        reaction = self._built(shape)
        reaction.entry = entry
        assert reaction.copy().entry is entry

    def test_the_two_copies_are_one_implementation(self):
        """
        `TemplateReaction.copy` and `LibraryReaction.copy` must not be two bodies that
        happen to agree -- that is how round 107's four sites drifted in the first place.
        Both are one call to the shared helper, differing only in the policy argument.

        **Structural**, but it runs at `56156ac9a` and fails on the base's actual source
        -- two bodies, each with its own `deepcopy` calls -- rather than on a missing
        name. It pins the repair's shape; the four `copy` arms above carry the defect.
        """
        for site in (TemplateReaction.copy, LibraryReaction.copy):
            body = inspect.getsource(site).split('"""')[-1]
            assert "copy_reaction(" in body, (
                "a copy() has stopped delegating to the shared helper")
            assert "deepcopy" not in body, (
                "a copy() deepens a field itself again; that is the per-field memo this "
                "round removed")

    def test_molecule_copy_still_preserves_atom_order(self):
        """
        The repair does not rely on this, but `Molecule.copy` does -- it pairs
        `self.vertices[i]` with `other.vertices[i]` to carry connectivity across. Pinned
        here because if that stops holding, a whole class of index-based reasoning about
        copied structures goes quietly wrong.

        **Green at `56156ac9a`** -- an upstream contract, pinned, not repaired.
        """
        molecule = Molecule(smiles="C[CH2]")
        other = molecule.copy(deep=True)
        assert [a.element.symbol for a in other.atoms] == \
               [a.element.symbol for a in molecule.atoms]
        assert all(a is not b for a in other.atoms for b in molecule.atoms)

    def test_a_deepcopy_cannot_do_this_and_the_reason_is_upstream(self):
        """
        The negative control on the choice of mechanism, so the workaround cannot outlive
        its reason.

        `copy()` reproduces the deepened half through `pickle` rather than `deepcopy`,
        which looks like the odd choice until the obvious one is tried:
        `Molecule.__deepcopy__` is ``return self.copy(deep=True)`` -- it accepts the memo
        and discards it, so no amount of memo-sharing by a caller can make two references
        to one molecule come out as two references to one copy. If upstream ever honours
        the memo, this fails and says the simpler mechanism has become available.

        **Green at `56156ac9a`** -- it measures upstream, which this round did not touch.
        """
        from copy import deepcopy

        molecule = Molecule(smiles="CC")
        atom = molecule.atoms[0]
        memo = {}
        copied_atom = deepcopy(atom, memo)
        copied_molecule = deepcopy(molecule, memo)
        assert not any(copied_atom is a for a in copied_molecule.atoms), (
            "Molecule.__deepcopy__ now honours the memo, so copy() could use a single "
            "deepcopy instead of a pickle round trip -- simplify it")


# ---------------------------------------------------------------------------------------
# Round 110
# ---------------------------------------------------------------------------------------


@pytest.mark.database
class TestTheTransportCarriesThePayloadToo:
    """
    Round 110's first HIGH. Round 108 chose the transport for one property -- one memo, so
    two fields referring to one object still refer to one object on the other side -- and
    it silently lost another. `pickle` delegates to each class's own ``__reduce__``, and
    ``Atom.__reduce__`` (``rmgpy/molecule/molecule.py:140``) never mentions ``id``,
    ``coords`` or ``props``. The body it replaced called `Species.copy(deep=True)`, which
    reproduces all three.

    Measured at `2e4ff991d` on a reaction this family generated, from molecules that had
    been through production's own `generate_resonance_structures()`: ids
    ``-32768…-32765`` came back ``-1`` and ``props {'inRing': False}`` came back ``{}``.
    Atom ids drive resonance-structure correspondence and ``'inRing'`` feeds group
    matching, so this is a wrong answer rather than a crash.

    **Round 108's tests cannot see it.** They compare `pairs` members and labelled atoms
    by identity -- correctly, that is their property -- and never look at what the objects
    those references point at actually hold. The two properties are asserted on one copy
    here for that reason.

    The fixture is the family's own output throughout. Nothing under test is assigned by
    the test.
    """

    @classmethod
    def setup_class(cls):
        families_path = os.path.join(settings["database.directory"], "kinetics", "families")
        if not os.path.isdir(os.path.join(families_path, REAL_FAMILY)):
            pytest.skip(f"database at {settings['database.directory']} has no "
                        f"{REAL_FAMILY} to generate reactions from")
        database = KineticsDatabase()
        database.load_families(path=families_path, families=[REAL_FAMILY])
        cls.family = database.families[REAL_FAMILY]

    def _generated(self, delete_labels=True):
        reactants = []
        for smiles in ("[Li+]", "[CH3]"):
            species = Species(molecule=[Molecule(smiles=smiles)])
            species.generate_resonance_structures()
            reactants.append(species.molecule[0])
        reactions = self.family.generate_reactions(reactants, delete_labels=delete_labels)
        assert reactions, "the family generated nothing, so there is nothing to copy"
        return reactions[0]

    @staticmethod
    def _atoms_of(reaction):
        atoms = []
        for structure in list(reaction.reactants) + list(reaction.products):
            for molecule in (structure.molecule if isinstance(structure, Species)
                             else [structure]):
                atoms.extend(molecule.atoms)
        return atoms

    def test_the_fixture_carries_the_state_before_anything_is_copied(self):
        """
        Anti-vacuity, and the answer to "did the test assign what it then checks".

        **Green at `2e4ff991d`**: it is the guard on the tests below, not evidence of the
        defect. If the family ever stops assigning atom ids this goes red and says so,
        rather than letting the assertions below pass over ``-1 == -1``.
        """
        atoms = self._atoms_of(self._generated(delete_labels=False))
        assert any(atom.id != -1 for atom in atoms), (
            "no atom the family generated carries an assigned id, so a test that the copy "
            "keeps ids would pass whether or not the copy keeps anything")
        assert any(atom.props for atom in atoms), (
            "no atom the family generated carries props, same problem")

    def test_the_copy_keeps_both_the_references_and_the_atom_state(self):
        """
        The acceptance for this round's HIGH 1, on **one** copy.

        Both properties or neither: the references *between* the fields (round 108) and
        the state inside the objects those references point at (round 110). Splitting them
        across two tests is what let a mechanism that satisfied one and broke the other
        land, so they are asserted together here.

        **Behavioural** at `2e4ff991d`: the ids and props assertions fail there, and the
        identity assertions pass there and must keep passing.
        """
        reaction = self._generated(delete_labels=False)
        before = self._atoms_of(reaction)
        copy = reaction.copy()
        after = self._atoms_of(copy)

        assert len(after) == len(before)
        for original, copied in zip(before, after):
            assert copied is not original, "the copy is aliasing the original's atoms"
            assert copied.id == original.id, (
                "atom id {0} came back as {1}; ids drive resonance-structure "
                "correspondence".format(original.id, copied.id))
            assert copied.props == original.props, (
                "atom props {0!r} came back as {1!r}; 'inRing' feeds group "
                "matching".format(original.props, copied.props))

        owned = list(copy.reactants) + list(copy.products)
        for pair in copy.pairs:
            for member in pair:
                assert any(member is s for s in owned), (
                    "round 108's property broke: a pairs member is a species the copy "
                    "does not own")
        for group in copy.labeled_atoms.values():
            for label, atom in group.items():
                for one in (atom if isinstance(atom, list) else [atom]):
                    assert any(one is a for a in after), (
                        "round 108's property broke: labelled atom {0!r} is inside no "
                        "molecule the copy owns".format(label))

    def test_the_familys_own_output_can_be_copied_at_all(self):
        """
        `generate_reactions()` **deletes** `labeled_atoms` once the labels have been read
        back (`family.py:2651`), so every reaction it returns is missing a field the copy
        policy lists. Round 108 applied the deepened set from the table rather than from
        the object, so `copy()` raised `ReactionStateNotCarried` on the family's own
        output -- the commonest reaction in the codebase.

        **Behavioural** at `2e4ff991d`: it raises there.
        """
        reaction = self._generated()
        assert not hasattr(reaction, "labeled_atoms"), (
            "the family no longer deletes labeled_atoms, so this test is no longer about "
            "the case it was written for")
        copy = reaction.copy()
        assert type(copy) is type(reaction)
        assert not hasattr(copy, "labeled_atoms"), (
            "the copy invented a field the original did not have")

    def test_the_copied_atoms_keep_the_interned_atom_types(self):
        """
        Why `Atom` keeps its own reducer and has what it drops added back, rather than
        being carried wholesale.

        `Atom.__reduce__` stores ``atomtype.label`` and restores ``ATOMTYPES[label]``, so
        a copied atom's type is the *interned* object. `AtomType` inherits identity
        equality and `is_specific_case_of` is a membership test over those objects, so a
        reducer that carried the type by value would trade three lost fields for silent
        group mismatches.

        **Green at `2e4ff991d`** -- it pins a property the repair had to preserve, not one
        it added.
        """
        reaction = self._generated(delete_labels=False)
        copy = reaction.copy()
        for original, copied in zip(self._atoms_of(reaction), self._atoms_of(copy)):
            assert copied.atomtype is original.atomtype, (
                "the copy's atom type is a different object from the interned one")

    def test_deepcopy_of_a_reaction_is_the_copy(self):
        """
        `deepcopy(reaction)` is a live path -- `family.py:3854` and
        `rmgpy/data/kinetics/database.py:755` both take one -- and it recurses into
        `Molecule.__deepcopy__`, which discards the memo, so the labelled atoms come back
        detached.

        **Behavioural** at `2e4ff991d`: the labelled atoms are severed there.
        """
        from copy import deepcopy

        reaction = self._generated(delete_labels=False)
        deep = deepcopy(reaction)
        assert type(deep) is type(reaction), (
            "deepcopy lost the subclass, which it did not do before")
        atoms = self._atoms_of(deep)
        for group in deep.labeled_atoms.values():
            for label, atom in group.items():
                for one in (atom if isinstance(atom, list) else [atom]):
                    assert any(one is a for a in atoms), (
                        "after deepcopy, labelled atom {0!r} is inside no molecule the "
                        "copy owns".format(label))


class TestTheClassesWhoseOwnReducerLoses:
    """
    The rest of round 110's first HIGH, and the whole of its second, on fixtures that do
    not need the database.

    `Fragment` and `CuttingLabel` inherit reducers that name `Molecule` and `Atom` as the
    class to rebuild; the cutting label's symbol is then handed to `get_element` and
    raises ``KeyError: 'R'``. `Molecule.__reduce__` passes ``metal`` and ``facet`` into
    ``__init__``'s ``inchi`` and ``smiles`` -- the fifth and sixth positional parameters
    rather than the seventh and eighth. `Species.__reduce__` omits four fields. All three
    were reached only because round 108 made `copy()` go through those reducers; the deep
    copy it replaced handled all three.
    """

    @staticmethod
    def _reaction(reactants, products):
        reaction = TemplateReaction(reactants=reactants, products=products,
                                    family="A_Family")
        reaction.labeled_atoms = {"reactants": {}, "products": {}}
        return reaction

    def test_a_fragment_reaction_can_be_copied(self):
        """
        Round 110's HIGH 2, behavioural: ``KeyError: 'R'`` at `2e4ff991d`.

        A working feature -- `deepcopy` of a fragment succeeds at both tips, which the
        negative control below pins -- broken by a change to how reactions are copied.
        """
        from rmgpy.molecule.fragment import CuttingLabel, Fragment

        reaction = self._reaction([Fragment().from_smiles_like_string("CCR")],
                                  [Fragment().from_smiles_like_string("[CH3]")])
        copy = reaction.copy()

        reactant = copy.reactants[0]
        assert isinstance(reactant, Fragment), (
            "the copied fragment came back as a {0}".format(type(reactant).__name__))
        labels = [atom for atom in reactant.atoms if isinstance(atom, CuttingLabel)]
        assert labels, "the copy has no CuttingLabel; the cutting label became an atom"
        assert [atom.name for atom in labels] == ["R"], (
            "the cutting label lost its name, which is the only thing that identifies it")
        assert all(atom is not other
                   for atom, other in zip(reactant.atoms,
                                          reaction.reactants[0].atoms)), (
            "the copy is aliasing the original's atoms")

    def test_the_deep_copy_this_replaced_could_do_it(self):
        """
        Anti-vacuity for the test above: fragment copying is a feature that worked, not
        one that never did.

        **Green at `2e4ff991d`** -- it measures the path round 108 replaced.
        """
        from copy import deepcopy

        from rmgpy.molecule.fragment import CuttingLabel, Fragment

        fragment = Fragment().from_smiles_like_string("CCR")
        deep = deepcopy(fragment)
        assert isinstance(deep, Fragment)
        assert any(isinstance(atom, CuttingLabel) for atom in deep.atoms)

    def test_a_surface_reaction_can_be_copied(self):
        """
        Behavioural: ``KeyError: 'Pt'`` at `2e4ff991d`, because ``metal='Pt'`` is read as
        an InChI. Every surface reaction in the codebase copies through this path.
        """
        molecule = Molecule(smiles="CC")
        molecule.metal, molecule.facet = "Pt", "111"
        reaction = self._reaction([molecule], [Molecule(smiles="[CH3]")])

        copy = reaction.copy()
        assert (copy.reactants[0].metal, copy.reactants[0].facet) == ("Pt", "111")

    def test_a_copied_species_keeps_what_its_own_reducer_omits(self):
        """
        Behavioural: `Species.__reduce__` passes ten fields and the class holds more.
        `symmetry_number` is the one that matters most -- it divides into every rate the
        species appears in -- and it came back ``-1.0`` from ``2.0``.
        """
        reactant = Species(label="ethane", molecule=[Molecule(smiles="CC")])
        reactant.symmetry_number = 2
        reactant.aug_inchi = "an augmented inchi"
        reactant.creation_iteration = 3
        reactant.explicitly_allowed = True
        reaction = self._reaction(
            [reactant], [Species(label="ethyl", molecule=[Molecule(smiles="C[CH2]")])])

        copy = reaction.copy().reactants[0]
        assert copy is not reactant
        assert copy.symmetry_number == 2
        assert copy.aug_inchi == "an augmented inchi"
        assert copy.creation_iteration == 3
        assert copy.explicitly_allowed is True

    def test_every_class_that_loses_state_is_registered(self):
        """
        The census, mechanised. Every class reachable in a reaction's deepened state is
        round-tripped through its *own* reducer and asked what it lost; anything that lost
        something must be in `_LOSSY_REDUCERS` with a reason. `Bond`, `TransitionState`
        and the kinetics models lose nothing and are deliberately absent.

        This is what keeps the registry from going stale in the direction that matters: a
        class that starts losing state, or one that stops.

        **Structural** at `2e4ff991d`, which has no `_LOSSY_REDUCERS`.
        """
        import numpy as np

        from rmgpy.data.kinetics.family import _LOSSY_REDUCERS, object_state
        from rmgpy.molecule.fragment import CuttingLabel, Fragment
        from rmgpy.molecule.molecule import Atom

        plain = (type(None), bool, int, float, complex, str, bytes, np.ndarray)

        def is_plain(value):
            if isinstance(value, plain):
                return True
            if isinstance(value, (list, tuple, set, frozenset)):
                return all(is_plain(item) for item in value)
            if isinstance(value, dict):
                return all(is_plain(k) and is_plain(v) for k, v in value.items())
            return False

        def loses(obj):
            try:
                back = pickle.loads(pickle.dumps(obj, pickle.HIGHEST_PROTOCOL))
            except Exception:                                    # noqa: BLE001
                return True
            if type(back) is not type(obj):
                return True
            for name, value in object_state(obj).items():
                if not is_plain(value):
                    continue
                other = getattr(back, name, None)
                if isinstance(value, np.ndarray) or isinstance(other, np.ndarray):
                    if not np.array_equal(value, other):
                        return True
                elif not (value is other or value == other):
                    return True
            return False

        molecule = Molecule(smiles="CC")
        atom = molecule.atoms[0]
        atom.id, atom.props = -424242, {"inRing": False}
        surface = Molecule(smiles="CC")
        surface.metal, surface.facet = "Pt", "111"
        species = Species(label="x", molecule=[Molecule(smiles="CC")])
        species.symmetry_number = 2
        fragment = Fragment().from_smiles_like_string("CCR")
        cutting = [a for a in fragment.atoms if isinstance(a, CuttingLabel)][0]

        reachable = {
            Atom: atom,
            Molecule: surface,
            Species: species,
            Fragment: fragment,
            CuttingLabel: cutting,
            type(list(atom.edges.values())[0]): list(atom.edges.values())[0],
            Arrhenius: Arrhenius(A=(1.0, "s^-1"), n=0, Ea=(0, "kJ/mol"), comment="c"),
            TransitionState: TransitionState(label="ts"),
        }
        for cls, instance in reachable.items():
            if loses(instance):
                assert cls in _LOSSY_REDUCERS, (
                    "{0} loses state through its own reducer and no complete reducer is "
                    "registered for it".format(cls.__name__))
            else:
                assert cls not in _LOSSY_REDUCERS, (
                    "{0} no longer loses anything through its own reducer, so the "
                    "registered workaround has outlived its reason".format(cls.__name__))
        for cls, (_, reason) in _LOSSY_REDUCERS.items():
            assert reason and len(reason) > 30, (
                "{0} is registered without a reason worth reading".format(cls.__name__))

    def test_the_probe_that_derives_the_dropped_set_can_see_every_field(self):
        """
        What makes the completeness of the `Atom` restore *testable* rather than asserted.

        The restore set is measured: a probe atom carrying a distinguishable value in
        every plain field is round-tripped, and whatever comes back changed is what gets
        restored. A field with no distinguishable value in the probe would be classified
        as surviving whether it survives or not -- so this asserts the probe actually
        distinguishes every plain field `Atom` has. An `Atom` that grows one turns this
        red and names it.

        **Structural** at `2e4ff991d`, which has no probe.
        """
        import numpy as np

        from rmgpy.data.kinetics.family import (_atom_probe, fields_the_reducer_drops,
                                                _is_plain, writable_fields)
        from rmgpy.molecule.molecule import Atom

        #: Fields the probe cannot give a distinguishable value to, each because its value
        #: is a reference rather than data: `edges` is a dict keyed by the atoms this atom
        #: is bonded to, and `atomtype` and `mapping` are objects that happen to be None
        #: on a bare atom. The derivation's answer for these is not evidence either way,
        #: and round 108's `is` assertions are what audits them. Anything *else* that
        #: cannot be distinguished is a hole in the derivation.
        references = {"edges", "atomtype", "mapping"}

        probe = _atom_probe()
        default = Atom(element="C")
        undistinguished = []
        for name in sorted(writable_fields(Atom) - references):
            mine, theirs = getattr(probe, name), getattr(default, name)
            if not _is_plain(mine):
                continue
            if isinstance(mine, np.ndarray) or isinstance(theirs, np.ndarray):
                if np.array_equal(mine, theirs):
                    undistinguished.append(name)
            elif mine == theirs:
                undistinguished.append(name)
        assert not undistinguished, (
            "these plain fields of Atom carry the default value in the probe, so the "
            "derivation cannot tell whether the reducer drops them: {0}".format(
                undistinguished))
        assert not references - writable_fields(Atom), (
            "these names are excused from the probe but are no longer writable fields of "
            "Atom, so the excuse has gone stale: {0}".format(
                sorted(references - writable_fields(Atom))))
        assert fields_the_reducer_drops(Atom) >= {"id", "coords", "props"}, (
            "the derivation no longer finds the three fields the review measured; either "
            "upstream fixed Atom.__reduce__ -- in which case say so and delete the "
            "workaround -- or the derivation has stopped working")


class TestThePartitionHasNoSilentDefault:
    """
    Round 110's first MEDIUM. `copy_reaction`'s two tables were complementary *for the
    fields somebody had classified*; for everything else the default was to carry by
    reference, and the default was silent. `TemplateReaction.template` is a mutable list
    and `specific_collider` is a `Species`, and neither was in either table, so
    ``reaction.copy().template is reaction.template`` was ``True``.

    There is no default now: a field in neither table is refused by name.
    """

    @staticmethod
    def _reaction(**fields):
        reaction = TemplateReaction(
            reactants=[Species(label="a", molecule=[Molecule(smiles="CC")])],
            products=[Species(label="b", molecule=[Molecule(smiles="C[CH2]")])],
            family="A_Family", **fields)
        reaction.labeled_atoms = {"reactants": {}, "products": {}}
        return reaction

    def test_the_template_is_not_shared_with_the_copy(self):
        """Behavioural: ``True`` at `2e4ff991d`, on a list the family rewrites in place."""
        reaction = self._reaction(template=["C/H3/Cs\\H3", "O_pri_rad"])
        copy = reaction.copy()
        assert copy.template == reaction.template
        assert copy.template is not reaction.template, (
            "the copy shares the original's template list, so appending to one appends "
            "to the other")

    def test_the_specific_collider_is_not_shared_with_the_copy(self):
        """Behavioural: the same, on a `Species` rather than a list."""
        collider = Species(label="Ar", molecule=[Molecule(smiles="[Ar]")])
        reaction = self._reaction(specific_collider=collider)
        copy = reaction.copy()
        assert copy.specific_collider is not None
        assert copy.specific_collider is not collider, (
            "the copy shares the original's collider species")
        assert copy.specific_collider.label == "Ar"

    def test_a_field_nobody_classified_is_refused_by_name(self):
        """
        The general form, and the reason this is a partition rather than a filter.

        **Behavioural** at `2e4ff991d`: no exception there, and the field is aliased.
        `ReactionStateNotCarried` is the base class of the exception the repair raises, so
        this test names only what the base already had and its red state is the missing
        refusal rather than a missing import.
        """
        from rmgpy.data.kinetics.family import ReactionStateNotCarried

        reaction = self._reaction()
        reaction.a_field_nobody_classified = ["a mutable one"]
        with pytest.raises(ReactionStateNotCarried) as raised:
            reaction.copy()
        assert "a_field_nobody_classified" in str(raised.value), (
            "the refusal does not name the field, so whoever hits it cannot act on it")

    def test_the_two_tables_are_a_partition_of_what_a_reaction_holds(self):
        """
        Every field of every shape lands in exactly one table, checked against the objects
        rather than against a list.

        **Structural** at `2e4ff991d`, which has no `_COPIED_BY_REFERENCE`.
        """
        from rmgpy.data.kinetics.family import (_COPIED_BY_REFERENCE,
                                                _NOT_COPIED_BY_REFERENCE,
                                                _TEMPLATE_NOT_COPIED_BY_REFERENCE,
                                                state_fields)

        shapes = ((TemplateReaction(reactants=[], products=[], family="F"),
                   _TEMPLATE_NOT_COPIED_BY_REFERENCE),
                  (LibraryReaction(reactants=[], products=[], library="L"),
                   _NOT_COPIED_BY_REFERENCE))
        for reaction, deepened in shapes:
            held = state_fields(reaction)
            unclassified = held - set(deepened) - set(_COPIED_BY_REFERENCE)
            assert not unclassified, (
                "{0} holds fields in neither table: {1}".format(
                    type(reaction).__name__, sorted(unclassified)))
            both = held & set(deepened) & set(_COPIED_BY_REFERENCE)
            assert not both, (
                "{0} holds fields in both tables: {1}".format(
                    type(reaction).__name__, sorted(both)))
        for name in ("template", "specific_collider"):
            assert (name in _TEMPLATE_NOT_COPIED_BY_REFERENCE
                    and name not in _COPIED_BY_REFERENCE), (
                "{0!r} is mutable and must be deepened, not aliased".format(name))
        for name, reason in _COPIED_BY_REFERENCE.items():
            assert reason, "{0!r} is shared without a reason".format(name)


class TestTheManifestIdentityBracketsTheRead:
    """
    Round 110's second MEDIUM. `os.fstat(fd)` came before `handle.read()` on the same
    descriptor. The descriptor defeats a *rename* between the two -- that is what closed
    round 99 -- and does nothing about a writer that rewrites the same file **in place**
    while the read is in flight: the identity is then the old version's and the content
    the new one's, and `resolve_quarantine` caches the second under the first.
    """

    MANIFEST = ("state = 'QUARANTINED FOR QUANTITATIVE PLASMA USE'\n"
                "reason = '{0}'\n"
                "kinetics_class = 'Marcus'\n")

    def _family(self, tmp_path):
        family_path = tmp_path / "kinetics" / "families" / "A_Family"
        family_path.mkdir(parents=True)
        manifest = family_path / QUARANTINE_FILENAME
        manifest.write_text(self.MANIFEST.format("the reason the first version gave"))
        return str(family_path), manifest

    def test_a_manifest_rewritten_between_the_two_reads_is_refused(self, tmp_path,
                                                                   monkeypatch):
        """
        The acceptance, and the only way to state it: rewrite the file *between* the
        identity call and the content call and require a refusal.

        **Behavioural** at `2e4ff991d`: the content comes back as the rewritten version
        keyed to the identity of the version before it.
        """
        from rmgpy.data.kinetics.quarantine import (QUARANTINE_UNREADABLE,
                                                    _read_manifest)

        family_path, manifest = self._family(tmp_path)
        rewritten = self.MANIFEST.format("a reason nobody has ever approved")
        real_fstat = os.fstat
        state = {"done": False}

        def rewriting_fstat(fd):
            info = real_fstat(fd)
            if not state["done"] and os.path.samestat(info, os.stat(str(manifest))):
                state["done"] = True
                manifest.write_text(rewritten)          # in place: same inode, new bytes
            return info

        monkeypatch.setattr(os, "fstat", rewriting_fstat)
        content, identity = _read_manifest(family_path)

        assert state["done"], (
            "the rewrite never fired, so this test did not exercise the window it was "
            "written for")
        assert content is QUARANTINE_UNREADABLE and identity is None, (
            "the manifest changed under the read and was answered for anyway: the "
            "content is {0!r}".format(content if content is QUARANTINE_UNREADABLE
                                      else (content or "")[:60]))
        assert content is not None, (
            "round 111: a refusal must not be spelled the same way as an absent "
            "manifest -- that is what made `family.quarantine = None` mean both 'this "
            "family is fine' and 'I could not check'")

    def test_an_unchanged_manifest_is_still_read(self, tmp_path):
        """
        Anti-vacuity: the refusal must cost nothing in the ordinary case, or it would be
        indistinguishable from breaking the loader.

        **Green at `2e4ff991d`**.
        """
        from rmgpy.data.kinetics.quarantine import _read_manifest

        family_path, manifest = self._family(tmp_path)
        content, identity = _read_manifest(family_path)
        assert content is not None and "the reason the first version gave" in content
        assert identity is not None

    def test_the_identity_is_taken_from_the_descriptor_on_both_sides(self):
        """
        The shape, so a later repair cannot satisfy the test above by re-stat'ing the
        *name* -- which would reintroduce the two-resolutions defect round 95 closed.

        **Structural** at `2e4ff991d`: one `os.fstat` call there, two here.
        """
        from rmgpy.data.kinetics import quarantine as module

        body = inspect.getsource(module._read_manifest)
        assert body.count("os.fstat(") == 2, (
            "the identity is taken {0} time(s) from the descriptor; it must be taken "
            "before and after the read".format(body.count("os.fstat(")))
        assert "os.stat(" not in body, (
            "the identity is being taken from a path again somewhere in _read_manifest; "
            "a name is resolved afresh every time it is used and that is the defect "
            "round 95 closed")


class TestTheFlagTheReviewNamed:
    """
    The addendum, routed from another ticket: `allow_max_rate_violation` is lost in
    production.

    Two sites were named. On this branch one of them is **already closed** and has been
    since round 102: `get_library_reactions` does omit the flag from its constructor calls
    and `_carry_entry_fields` supplies it on the next line, from a field set discovered
    from `Reaction` rather than written down. The other is live and is fixed here.
    """

    def test_the_loader_carries_the_flag_despite_the_constructor_omitting_it(self):
        """
        **Green at `2e4ff991d`**, and recorded rather than claimed: the review is right
        about the constructor call and right for any branch without round 102's carry.
        """
        item = Reaction(reactants=[Species(label="a", molecule=[Molecule(smiles="CC")])],
                        products=[Species(label="b", molecule=[Molecule(smiles="C[CH2]")])])
        item.allow_max_rate_violation = True
        entry = Entry(index=1, label="a <=> b", item=item,
                      data=Arrhenius(A=(1.0, "s^-1"), n=0, Ea=(0, "kJ/mol")))
        library = KineticsLibrary(label="L")
        library.entries = {"a <=> b": entry}
        library.auto_generated = False

        reaction = library.get_library_reactions()[0]
        assert reaction.allow_max_rate_violation is True
        source = inspect.getsource(KineticsLibrary.get_library_reactions)
        assert "allow_max_rate_violation" not in source, (
            "the constructor calls now name the flag, so this test's premise -- that the "
            "derived carry is what supplies it -- has changed")

    def test_the_base_reaction_copy_carries_the_flag(self):
        """Behavioural: ``False`` at `2e4ff991d`, and the default is the value that hides
        the loss."""
        reaction = Reaction(
            reactants=[Species(label="a", molecule=[Molecule(smiles="CC")])],
            products=[Species(label="b", molecule=[Molecule(smiles="C[CH2]")])])
        reaction.allow_max_rate_violation = True
        reaction.rank = 5
        reaction.is_forward = True

        copy = reaction.copy()
        assert copy.allow_max_rate_violation is True
        assert copy.rank == 5, "`rank` was beside it and missing the same way"
        assert copy.is_forward is True, "`is_forward` too"
        assert copy.k_effective_cache == {}, (
            "__new__ leaves a cdef public dict unset and reading one raises, so a copy "
            "without this cannot be asked for a rate coefficient")

    def test_the_base_reaction_reduce_carries_the_flag(self):
        """Behavioural: a deepcopy goes through `__reduce__`, and it dropped the flag."""
        from copy import deepcopy

        reaction = Reaction(
            reactants=[Species(label="a", molecule=[Molecule(smiles="CC")])],
            products=[Species(label="b", molecule=[Molecule(smiles="C[CH2]")])])
        reaction.allow_max_rate_violation = True
        reaction.is_forward = False

        for after in (pickle.loads(pickle.dumps(reaction)), deepcopy(reaction)):
            assert after.allow_max_rate_violation is True
            assert after.is_forward is False

    def test_the_base_reaction_copy_keeps_its_pairs(self):
        """
        Behavioural, and the reason `__deepcopy__` could be added at all: while
        `Reaction.copy` severed `pairs`, routing `deepcopy` through `copy()` would have
        been a regression rather than a repair.
        """
        reaction = Reaction(
            reactants=[Species(label="a", molecule=[Molecule(smiles="CC")])],
            products=[Species(label="b", molecule=[Molecule(smiles="C[CH2]")])])
        reaction.pairs = [(reaction.reactants[0], reaction.products[0])]

        copy = reaction.copy()
        owned = list(copy.reactants) + list(copy.products)
        for pair in copy.pairs:
            for member in pair:
                assert any(member is s for s in owned), (
                    "a pairs member is a species the copy does not own; "
                    "Species.__eq__ is identity, so reactants.index(pair[0]) raises")
                assert not any(member is s for s in
                               list(reaction.reactants) + list(reaction.products)), (
                    "the copy's pairs still point at the original's species")

    def test_the_flag_is_named_in_the_partition_rather_than_carried_by_default(self):
        """
        What the addendum asked for beyond the two fixes: the flag is classified
        explicitly, so the next person reading the table sees it.

        **Structural** at `2e4ff991d`.
        """
        from rmgpy.data.kinetics.family import _COPIED_BY_REFERENCE

        assert "allow_max_rate_violation" in _COPIED_BY_REFERENCE
        assert "Reaction.copy" in _COPIED_BY_REFERENCE["allow_max_rate_violation"], (
            "the entry does not say why it is worth naming")
        assert "is_forward" in _COPIED_BY_REFERENCE


IONISATION_FAMILY = "Plasma_Electron_Impact_Ionization"


def _payload_of(reaction):
    """
    Every atom's ``(id, props)`` in a reaction, in order.

    The payload, not the labels. Round 110's pickle assertions collapse structures to
    labels (`TestEveryTransformReproducesTheWholeState`), which is a comparison no severed
    payload can fail: a reaction whose atoms have all lost their ids still prints the same
    SMILES.
    """
    payload = []
    for structure in list(reaction.reactants) + list(reaction.products):
        molecules = (structure.molecule if isinstance(structure, Species)
                     else [structure])
        for molecule in molecules:
            for atom in molecule.atoms:
                payload.append((atom.id, dict(atom.props)))
    return payload


@pytest.mark.database
class TestTheOtherTransportGetsTheSameReducers:
    """
    Round 111's first HIGH. Round 110 completed the five lossy reducers and installed
    them **inside `copy()`**, as one pickler's ``dispatch_table``. A dispatch table
    belongs to the pickler carrying it, so every other transport went on getting the
    lossy ones -- and the one production takes on every parallel run is `Pool.map`
    (``rmgpy/rmg/react.py``), which serialises with `multiprocessing`'s `ForkingPickler`
    in both directions.

    Measured at `13e3227b2`: a model generated with ``procnum=2`` came back with every
    atom id ``-1`` and every ``props`` empty, against the same generation at ``procnum=1``
    -- and a fragment or surface reaction raised ``KeyError`` on the parallel path where
    the serial path succeeded.

    Every reaction here comes out of the family, and `generate_reactions` is called with
    **default** arguments. Round 110's own verification passed ``delete_labels=False``,
    which made the state observable and the path unreachable.
    """

    @classmethod
    def setup_class(cls):
        families_path = os.path.join(settings["database.directory"], "kinetics", "families")
        if not os.path.isdir(os.path.join(families_path, REAL_FAMILY)):
            pytest.skip(f"database at {settings['database.directory']} has no "
                        f"{REAL_FAMILY} to generate reactions from")
        cls.database = KineticsDatabase()
        cls.database.load_families(path=families_path, families=[REAL_FAMILY])
        cls.family = cls.database.families[REAL_FAMILY]

    @staticmethod
    def _reactant_species():
        species = []
        for smiles in ("[Li+]", "[CH3]"):
            one = Species(molecule=[Molecule(smiles=smiles)])
            one.generate_resonance_structures()
            species.append(one)
        return species

    def _generated(self):
        """The family's own output, with **default** arguments. No keyword is passed."""
        reactants = [one.molecule[0] for one in self._reactant_species()]
        reactions = self.family.generate_reactions(reactants)
        assert reactions, "the family generated nothing, so there is nothing to transport"
        return reactions[0]

    @staticmethod
    def _through_multiprocessing(value):
        from multiprocessing.reduction import ForkingPickler
        return ForkingPickler.loads(ForkingPickler.dumps(value))

    def test_the_fixture_carries_a_payload_before_anything_is_transported(self):
        """Anti-vacuity: the ids and props under test are production's, not the test's."""
        payload = _payload_of(self._generated())
        assert payload, "no atoms were found to compare"
        assert any(identifier != -1 for identifier, _ in payload), (
            "no atom carries an id, so an assertion that ids survive could not fail")
        assert any(props for _, props in payload), (
            "no atom carries props, so an assertion that props survive could not fail")

    def test_the_transport_multiprocessing_uses_carries_the_payload(self):
        """
        Behavioural at `13e3227b2`: ids ``-32768...`` came back ``-1`` and
        ``props {'inRing': False}`` came back ``{}``.
        """
        reaction = self._generated()
        before = _payload_of(reaction)
        after = _payload_of(self._through_multiprocessing(reaction))
        assert after == before, (
            "the transport `Pool.map` uses dropped per-atom state that `copy()` keeps; "
            "ids {0} -> {1}".format([b[0] for b in before[:4]],
                                    [a[0] for a in after[:4]]))

    def test_a_parallel_generation_reproduces_the_serial_one(self):
        """
        The finding in production's own terms: the same model, generated serially and in
        parallel, compared on atom payload rather than on labels.

        Behavioural at `13e3227b2`. `react()` is the real entry point -- this starts a
        real `Pool` -- with the loaded database installed where `get_db` reads it.
        """
        import rmgpy.data.rmg
        from rmgpy.rmg.react import react

        class _Stub:
            pass

        stub = _Stub()
        stub.kinetics = self.database
        previous = getattr(rmgpy.data.rmg, "database", None)
        rmgpy.data.rmg.database = stub
        try:
            reactants = self._reactant_species()
            task = [((reactants[0], reactants[1]), [REAL_FAMILY])]
            serial = react(list(task), procnum=1)[0]
            parallel = react(list(task), procnum=2)[0]
        finally:
            rmgpy.data.rmg.database = previous

        assert serial and parallel, "the generation produced nothing to compare"
        assert len(serial) == len(parallel)
        assert _payload_of(parallel[0]) == _payload_of(serial[0]), (
            "a model generated in parallel differs from the same model generated "
            "serially: the atoms came back without the ids that drive "
            "resonance-structure correspondence and without the props that feed group "
            "matching")

    def test_a_fragment_reaction_survives_the_parallel_transport(self):
        """
        Behavioural at `13e3227b2`: ``KeyError: 'R'``. `CuttingLabel` inherits
        `Atom.__reduce__`, which hands its symbol to `get_element`.
        """
        from rmgpy.molecule.fragment import CuttingLabel, Fragment

        reaction = Reaction(reactants=[Fragment().from_smiles_like_string("CCR")],
                            products=[Fragment().from_smiles_like_string("[CH3]")])
        after = self._through_multiprocessing(reaction)
        assert isinstance(after.reactants[0], Fragment)
        assert any(isinstance(atom, CuttingLabel) for atom in after.reactants[0].atoms), (
            "the cutting label did not survive the transport")

    def test_a_surface_reaction_survives_the_parallel_transport(self):
        """
        Behavioural at `13e3227b2`: ``KeyError: 'Pt'``. `Molecule.__reduce__` passes
        ``metal`` into ``__init__``'s ``inchi`` parameter.
        """
        molecule = Molecule(smiles="C")
        molecule.metal = "Pt"
        molecule.facet = "111"
        after = self._through_multiprocessing(
            Reaction(reactants=[molecule], products=[Molecule(smiles="C")]))
        assert (after.reactants[0].metal, after.reactants[0].facet) == ("Pt", "111")

    def test_both_transports_carry_one_table_rather_than_two(self):
        """
        The anti-drift property, and the reason this is a repair to the mechanism rather
        than a second copy of it: `multiprocessing`'s pickler is handed the *same* dict
        `copy()`'s pickler carries, so a class added to `_LOSSY_REDUCERS` is carried by
        both without anyone remembering to do it twice.

        **Structural** at `13e3227b2` -- it names something the repair adds.
        """
        from multiprocessing.reduction import ForkingPickler

        from rmgpy.data.kinetics.family import (COMPLETE_REDUCERS, _CompletePickler,
                                                _LOSSY_REDUCERS)

        assert set(COMPLETE_REDUCERS) == set(_LOSSY_REDUCERS)
        assert _CompletePickler.dispatch_table is COMPLETE_REDUCERS
        for cls, reducer in COMPLETE_REDUCERS.items():
            assert ForkingPickler._extra_reducers.get(cls) is reducer, (
                "{0} is completed for copy() and not for the parallel "
                "path".format(cls.__name__))

    def test_the_module_that_needs_the_registration_makes_it(self):
        """
        Where the repair lives. `react.py` is the module whose correctness depends on the
        registration, so it calls for it rather than inheriting it from an import graph
        that may be reordered.

        **Structural** at `13e3227b2`.
        """
        import rmgpy.rmg.react

        source = inspect.getsource(rmgpy.rmg.react)
        assert "install_complete_reducers()" in source, (
            "react.py does not ask for the reducers it depends on")


@pytest.mark.database
class TestASpeciesOnBothSidesKeepsItsSide:
    """
    Round 111's second HIGH, on the campaign's own reaction.

    `Reaction.copy()` built one ``id(original) -> copy`` map over reactants and products
    together, so a `Species` on both sides had its reactant entry overwritten by the
    product one. Every electron-impact reaction is that shape:
    `electron_placement.py:_place_declared_electrons` appends the *same* canonical
    electron object to both sides, so ``Li + e- => Li+ + e- + e-`` holds one object three
    times.

    Measured at `13e3227b2` on the resolver's own view object: pair 0's reactant-side
    member was in ``copy.products`` and not in ``copy.reactants``, and
    ``reactants.index(pair[0])`` raised ``ValueError``.

    The existing test (`test_the_base_reaction_copy_keeps_its_pairs`) uses distinct
    species only, which is why it stayed green through all of this.
    """

    @classmethod
    def setup_class(cls):
        families_path = os.path.join(settings["database.directory"], "kinetics", "families")
        if not os.path.isdir(os.path.join(families_path, IONISATION_FAMILY)):
            pytest.skip(f"database at {settings['database.directory']} has no "
                        f"{IONISATION_FAMILY} to generate the electron-impact case from")
        database = KineticsDatabase()
        database.load_families(path=families_path, families=[IONISATION_FAMILY])
        cls.family = database.families[IONISATION_FAMILY]

    def _view(self):
        from rmgpy.electron_placement import resolve_electron_placement

        reaction = self.family.generate_reactions([Molecule(smiles="[Li]")])[0]
        reaction.kinetics = self.family.get_kinetics(
            reaction, template_labels=reaction.template,
            degeneracy=reaction.degeneracy)[0][0]
        reaction.ensure_species()
        electron = Species(label="e", molecule=[Molecule(smiles="e")])
        view = resolve_electron_placement(
            reaction, [electron] + list(reaction.reactants) + list(reaction.products))
        view.generate_pairs()
        return view

    def test_the_fixture_holds_one_species_on_both_sides(self):
        """
        Anti-vacuity. Without this, every assertion below passes on a reaction whose
        species happen to be distinct -- which is exactly how the defect survived.
        """
        view = self._view()
        shared = [r for r in view.reactants for p in view.products if r is p]
        assert shared, (
            "the fixture does not hold one Species object on both sides, so it cannot "
            "exercise the defect")
        assert view.pairs, "the fixture has no pairs to sever"

    def test_a_reactant_side_pair_member_is_one_of_the_copys_reactants(self):
        """Behavioural at `13e3227b2`: pair 0's member was in products only."""
        view = self._view()
        copied = view.copy()
        for index, pair in enumerate(copied.pairs):
            assert any(species is pair[0] for species in copied.reactants), (
                "pair {0}'s reactant-side member ({1!s}) is not one of the copy's "
                "reactants".format(index, pair[0]))
            assert any(species is pair[1] for species in copied.products), (
                "pair {0}'s product-side member ({1!s}) is not one of the copy's "
                "products".format(index, pair[1]))

    def test_indexing_a_pair_member_does_not_raise(self):
        """
        Behavioural at `13e3227b2`: ``ValueError``. This is how the defect actually
        surfaces -- `Species.__eq__` is identity, so the exception lands in unrelated
        code far from the copy that caused it.
        """
        copied = self._view().copy()
        for pair in copied.pairs:
            copied.reactants.index(pair[0])
            copied.products.index(pair[1])

    def test_the_copy_owns_none_of_the_originals_species(self):
        """The property the side-aware map must not have cost: it is still a deep copy."""
        view = self._view()
        copied = view.copy()
        originals = list(view.reactants) + list(view.products)
        for species in list(copied.reactants) + list(copied.products):
            assert not any(species is other for other in originals)

    def test_the_collider_is_deepened_rather_than_aliased(self):
        """
        Behavioural at `13e3227b2`: ``copy().specific_collider is
        reaction.specific_collider``. The docstring says "deep copy", and an edit to the
        copy's collider reached the original's.
        """
        reaction = self._view().copy()
        reaction.specific_collider = Species(label="M", molecule=[Molecule(smiles="[He]")])
        copied = reaction.copy()
        assert copied.specific_collider is not reaction.specific_collider
        assert copied.specific_collider.label == reaction.specific_collider.label

    def test_a_collider_that_is_a_participant_stays_one_object_with_it(self):
        """
        The other half of the same property: deepening the collider independently would
        hand the copy a third-body species that is not the reactant it is supposed to be,
        which is the `pairs` defect one field over.
        """
        view = self._view()
        reaction = view.copy()
        reaction.specific_collider = reaction.reactants[0]
        copied = reaction.copy()
        assert copied.specific_collider is copied.reactants[0]


class TestARefusedManifestIsNotAnAbsentOne:
    """
    Round 111's first MEDIUM, and round 99's finding one layer down.

    `_read_manifest` answered ``(None, None)`` for every refusal -- a link, a FIFO, a
    permission error, a read or `fstat` that failed, a file rewritten mid-read -- and
    `load_family_quarantine` handed that to ``family.quarantine``, where it is the same
    value every ordinary family carries. `check_quarantine(..., family=family)` then
    called it **answered**, so a check that failed produced a clean bill of health.
    """

    @staticmethod
    def _families(tmp_path):
        root = str(tmp_path)
        readable = os.path.join(root, "kinetics", "families", "Readable")
        unreadable = os.path.join(root, "kinetics", "families", "Unreadable")
        absent = os.path.join(root, "kinetics", "families", "Absent")
        for path in (readable, unreadable, absent):
            os.makedirs(path)
        write_manifest(readable)
        # A link is refused by the O_NOFOLLOW descent: the manifest is here and cannot be
        # read, which is the state that must not read as absence.
        target = os.path.join(root, "elsewhere.py")
        with open(target, "w") as handle:
            handle.write(MANIFEST)
        os.symlink(target, os.path.join(unreadable, QUARANTINE_FILENAME))
        return root, readable, unreadable, absent

    def test_the_three_answers_are_three_values(self, tmp_path):
        """
        Behavioural at `13e3227b2` for the middle one: it was ``None``, the same value as
        the right-hand one.
        """
        from rmgpy.data.kinetics.quarantine import QUARANTINE_UNREADABLE

        _root, readable, unreadable, absent = self._families(tmp_path)

        assert isinstance(load_family_quarantine("Readable", readable), KineticsQuarantine)
        assert load_family_quarantine("Unreadable", unreadable) is QUARANTINE_UNREADABLE
        assert load_family_quarantine("Absent", absent) is None

    def test_a_refused_manifest_is_unanswered_at_the_gate(self, tmp_path, monkeypatch):
        """
        Behavioural at `13e3227b2`: `check_quarantine` returned in silence, which is what
        it does for a family that is genuinely not quarantined.
        """
        import rmgpy.data.rmg

        root, _readable, unreadable, _absent = self._families(tmp_path)
        _clear_gate_caches()
        monkeypatch.setitem(settings, "database.directory", root)
        monkeypatch.setattr(rmgpy.data.rmg, "database", None, raising=False)

        class _Family:
            label = "Unreadable"
            quarantine = load_family_quarantine("Unreadable", unreadable)

        reaction = make_library_reaction()
        reaction.kinetics = make_marcus()

        messages = []
        monkeypatch.setattr(logging, "warning",
                            lambda msg, *args, **kw: messages.append(msg % args))
        check_quarantine(reaction, "test", family=_Family())
        _clear_gate_caches()

        assert any("Cannot tell" in message for message in messages), (
            "a manifest that could not be read was answered for in silence")
        assert any("could not be read" in message for message in messages), (
            "the warning does not say which of the two unanswered cases this is; "
            '"load that family" is useless advice for a manifest that is already here')

    def test_an_absent_manifest_is_still_a_clean_answer(self, tmp_path, monkeypatch):
        """
        The control. The repair must not turn every ordinary family into an unanswered
        question -- that would put a warning into every run in the world.
        """
        import rmgpy.data.rmg

        root, _readable, _unreadable, absent = self._families(tmp_path)
        _clear_gate_caches()
        monkeypatch.setitem(settings, "database.directory", root)
        monkeypatch.setattr(rmgpy.data.rmg, "database", None, raising=False)

        class _Family:
            label = "Absent"
            quarantine = load_family_quarantine("Absent", absent)

        reaction = make_library_reaction()
        reaction.kinetics = make_marcus()

        messages = []
        monkeypatch.setattr(logging, "warning",
                            lambda msg, *args, **kw: messages.append(msg % args))
        check_quarantine(reaction, "test", family=_Family())
        _clear_gate_caches()

        assert not messages, "an ordinary family was reported as unanswered: {0}".format(
            messages[:1])

    def test_the_refusal_is_not_cached_as_an_answer(self, tmp_path, monkeypatch):
        """
        The condition that produces a refusal is transient. Caching it would make a
        momentary failure permanent for the life of the run -- and caching it as an
        *answer* would be the original defect with a memory.
        """
        import rmgpy.data.rmg

        root, _readable, unreadable, _absent = self._families(tmp_path)
        _clear_gate_caches()
        monkeypatch.setitem(settings, "database.directory", root)
        monkeypatch.setattr(rmgpy.data.rmg, "database", None, raising=False)

        assert resolve_quarantine("Unreadable") == (None, False)

        # Repair the manifest in place; the next lookup must see it rather than a cached
        # refusal.
        os.unlink(os.path.join(unreadable, QUARANTINE_FILENAME))
        write_manifest(unreadable)
        quarantine, answered = resolve_quarantine("Unreadable")
        _clear_gate_caches()

        assert answered is True
        assert isinstance(quarantine, KineticsQuarantine)


class TestTheEnumerationReadsTheDiskToo:
    """
    Round 111's third HIGH. `_warn_unattributable` enumerated **loaded** families only,
    so a library rate whose ``family:`` line is missing met an empty enumeration and
    passed with neither a refusal nor a warning -- while the quarantined family sat in the
    database directory, unloaded.

    That is the ordinary case rather than a contrived one:
    ``CoreEdgeReactionModel.add_seed_mechanism_to_core`` converts a seed reaction whose
    family is unavailable into a library reaction rather than loading the family, so
    "authorship lost, family not loaded" is the normal path for any foreign seed.
    """

    @staticmethod
    def _database_with_a_quarantine_on_disk(tmp_path, label="A_Family_On_Disk"):
        family_path = os.path.join(str(tmp_path), "kinetics", "families", label)
        os.makedirs(family_path)
        write_manifest(family_path)
        return str(tmp_path)

    @staticmethod
    def _unattributed_reaction():
        reaction = make_library_reaction(library="a_seed_with_no_family_line")
        reaction.kinetics = make_marcus()
        reaction.entry = Entry(index=1, label="no authorship", long_desc="nothing here")
        return reaction

    def test_an_unattributed_rate_is_reported_when_the_family_is_only_on_disk(
            self, tmp_path, monkeypatch):
        """
        Behavioural at `13e3227b2`: nothing was logged at all.
        """
        import rmgpy.data.rmg

        root = self._database_with_a_quarantine_on_disk(tmp_path)
        _clear_gate_caches()
        monkeypatch.setitem(settings, "database.directory", root)
        monkeypatch.setattr(rmgpy.data.rmg, "database", None, raising=False)

        messages = []
        monkeypatch.setattr(logging, "warning",
                            lambda msg, *args, **kw: messages.append(msg % args))
        check_quarantine(self._unattributed_reaction(), "test")
        _clear_gate_caches()

        assert messages, (
            "the rate was admitted in silence: the manifest is on disk, its criterion is "
            "this rate's own kinetics class, and the entry records no author")
        assert "A_Family_On_Disk" in messages[0], (
            "the warning does not name the family whose criterion matched")
        assert "NOT a refusal" in messages[0], (
            "the warning stopped saying that it is deliberately not a refusal")

    def test_nothing_is_said_when_the_database_quarantines_nothing(
            self, tmp_path, monkeypatch):
        """
        The control, and the promise this whole module rests on: a database with no
        manifest anywhere behaves exactly as it did before.
        """
        import rmgpy.data.rmg

        os.makedirs(os.path.join(str(tmp_path), "kinetics", "families", "Ordinary"))
        _clear_gate_caches()
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        monkeypatch.setattr(rmgpy.data.rmg, "database", None, raising=False)

        messages = []
        monkeypatch.setattr(logging, "warning",
                            lambda msg, *args, **kw: messages.append(msg % args))
        check_quarantine(self._unattributed_reaction(), "test")
        _clear_gate_caches()

        assert not messages, "an ordinary database gained a warning: {0}".format(
            messages[:1])

    def test_a_loaded_family_is_not_enumerated_twice(self, tmp_path, monkeypatch):
        """
        The loaded half and the disk half are the same families. Yielding one twice would
        double every count taken over this enumeration.
        """
        import rmgpy.data.rmg
        from rmgpy.data.kinetics.quarantine import iter_quarantines

        label = "A_Family_On_Disk"
        root = self._database_with_a_quarantine_on_disk(tmp_path, label)
        _clear_gate_caches()
        monkeypatch.setitem(settings, "database.directory", root)

        class _Family:
            def __init__(self, q):
                self.quarantine = q

        class _Kinetics:
            def __init__(self, f):
                self.families = f

        class _Database:
            def __init__(self, f):
                self.kinetics = _Kinetics(f)

        loaded = load_family_quarantine(
            label, os.path.join(root, "kinetics", "families", label))
        monkeypatch.setattr(rmgpy.data.rmg, "database",
                            _Database({label: _Family(loaded)}), raising=False)

        found = list(iter_quarantines())
        _clear_gate_caches()

        assert len(found) == 1, "the family was enumerated {0} times".format(len(found))
        # Not `found[0] is loaded`: that pinned round 111's stale path, where the cached
        # object was served instead of the manifest as it now reads (round 112's HIGH 4).
        assert (found[0].family_label, found[0].kinetics_class_name) == (label, "Marcus")

    def test_an_unreadable_manifest_is_not_enumerated_as_a_criterion(
            self, tmp_path, monkeypatch):
        """
        `QUARANTINE_UNREADABLE` has no ``applies_to``, so an enumeration that yielded it
        would raise inside the gate. It is not yielded; since round 112 the enumeration
        then raises `QuarantineEnumerationError`, because a manifest it could not read is a
        hole in it, not an absence.
        """
        import rmgpy.data.rmg
        from rmgpy.data.kinetics.quarantine import (QUARANTINE_UNREADABLE,
                                                    iter_quarantines)

        _clear_gate_caches()
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))

        class _Family:
            quarantine = QUARANTINE_UNREADABLE

        class _Kinetics:
            families = {"Unreadable": _Family()}

        class _Database:
            kinetics = _Kinetics()

        monkeypatch.setattr(rmgpy.data.rmg, "database", _Database(), raising=False)
        from rmgpy.data.kinetics.quarantine import QuarantineEnumerationError
        found = []
        with pytest.raises(QuarantineEnumerationError):
            # Round 112: an unreadable manifest might name any class, so the enumeration
            # is incomplete and says so, rather than ending as if there were nothing.
            for quarantine in iter_quarantines():
                found.append(quarantine)
        _clear_gate_caches()

        assert found == [], "a refusal was enumerated as if it were a criterion"


class TestTheEnumerationIsNeitherStaleNorSilent:
    """
    Round 112's fourth and fifth HIGH, both in `iter_quarantines`.

    HIGH 4: the loaded half yielded the cached family object and the disk half then
    skipped that family, so a manifest rewritten after load was enumerated as it used to
    read. HIGH 5: a families directory that could not be listed returned an empty
    enumeration, which `_warn_unattributable` could not tell from "nothing is quarantined".
    """

    ARRHENIUS_MANIFEST = MANIFEST.replace('"Marcus"', '"Arrhenius"')

    @staticmethod
    def _loaded_database(tmp_path, monkeypatch, label="F"):
        """A database directory holding one Marcus manifest, with that family loaded."""
        import rmgpy.data.rmg

        family_path = os.path.join(str(tmp_path), "kinetics", "families", label)
        os.makedirs(family_path)
        write_manifest(family_path)
        _clear_gate_caches()
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        loaded = load_family_quarantine(label, family_path)
        assert loaded.kinetics_class_name == "Marcus"
        _register_families(monkeypatch, {label: loaded})
        return family_path

    @staticmethod
    def _unattributed(kinetics):
        reaction = make_library_reaction(library="a_seed_with_no_family_line")
        reaction.kinetics = kinetics
        reaction.entry = Entry(index=1, label="no authorship", long_desc="nothing here")
        return reaction

    @staticmethod
    def _capture_warnings(monkeypatch):
        messages = []
        monkeypatch.setattr(logging, "warning",
                            lambda msg, *args, **kw: messages.append(msg % args))
        return messages

    @staticmethod
    def _unlistable(monkeypatch, path):
        """`path` denied as a mode-000 directory is, by patch; see `_deny_access`."""
        _deny_access(monkeypatch, path)
        with pytest.raises(PermissionError):
            os.listdir(path)

    def test_a_manifest_rewritten_after_load_is_enumerated_as_it_now_reads(
            self, tmp_path, monkeypatch):
        """
        The brief's sequence: load a Marcus manifest, rewrite it to Arrhenius on disk,
        enumerate. Behavioural at `71ae97bd5`: still ``[('F', 'Marcus')]``.
        """
        from rmgpy.data.kinetics.quarantine import iter_quarantines

        family_path = self._loaded_database(tmp_path, monkeypatch)
        write_manifest(family_path, self.ARRHENIUS_MANIFEST)

        found = [(q.family_label, q.kinetics_class_name) for q in iter_quarantines()]
        _clear_gate_caches()

        assert found == [("F", "Arrhenius")], (
            "the enumeration served the manifest as it read at load time: {0}".format(found))

    def test_an_unattributed_rate_of_the_rewritten_class_is_reported(
            self, tmp_path, monkeypatch):
        """
        What the stale enumeration cost: an unattributed Arrhenius rate admitted in silence
        while the manifest on disk names Arrhenius.
        """
        family_path = self._loaded_database(tmp_path, monkeypatch)
        write_manifest(family_path, self.ARRHENIUS_MANIFEST)
        messages = self._capture_warnings(monkeypatch)

        check_quarantine(self._unattributed(Arrhenius(A=(1e10, "m^3/(mol*s)"), n=0,
                                                      Ea=(10, "kJ/mol"))), "test")
        _clear_gate_caches()

        assert any("Arrhenius" in m and "family F (" in m for m in messages), (
            "the rewritten criterion was not applied: {0}".format(messages))

    def test_an_unlistable_families_directory_is_not_zero_quarantines(
            self, tmp_path, monkeypatch):
        """
        The brief's sequence: an unattributed library rate, an inaccessible families
        directory, then `check_quarantine`. Behavioural at `71ae97bd5`: no refusal and no
        warning, because the enumeration swallowed the `listdir` failure and came back
        empty.
        """
        import rmgpy.data.rmg

        families_root = os.path.join(str(tmp_path), "kinetics", "families")
        os.makedirs(os.path.join(families_root, "A_Family_On_Disk"))
        write_manifest(os.path.join(families_root, "A_Family_On_Disk"))
        _clear_gate_caches()
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        monkeypatch.setattr(rmgpy.data.rmg, "database", None, raising=False)
        messages = self._capture_warnings(monkeypatch)

        self._unlistable(monkeypatch, families_root)
        try:
            check_quarantine(self._unattributed(make_marcus()), "test")
        finally:
            _clear_gate_caches()

        assert messages, (
            "the rate was admitted in silence: the families directory could not be listed, "
            "and that was reported as there being no quarantine anywhere")
        assert any("could not" in m.lower() for m in messages), (
            "the warning does not say the enumeration was incomplete: {0}".format(messages))

    def test_an_unexaminable_family_directory_is_not_zero_quarantines(
            self, tmp_path, monkeypatch):
        """
        One layer down: the families directory lists, but one family's directory cannot
        be looked into. `os.path.lexists` answers False on that, which read as "this family
        carries no manifest".
        """
        import rmgpy.data.rmg

        family_path = os.path.join(str(tmp_path), "kinetics", "families", "A_Family_On_Disk")
        os.makedirs(family_path)
        write_manifest(family_path)
        _clear_gate_caches()
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        monkeypatch.setattr(rmgpy.data.rmg, "database", None, raising=False)
        messages = self._capture_warnings(monkeypatch)

        self._unlistable(monkeypatch, family_path)
        try:
            check_quarantine(self._unattributed(make_marcus()), "test")
        finally:
            _clear_gate_caches()

        assert any("could not" in m.lower() for m in messages), (
            "an unexaminable family read as an unquarantined one: {0}".format(messages))

    def test_the_public_enumeration_raises_rather_than_ending_early(
            self, tmp_path, monkeypatch):
        """An iterator that ends is a claim that there is nothing more; this one refuses."""
        import rmgpy.data.rmg
        from rmgpy.data.kinetics import quarantine as module

        families_root = os.path.join(str(tmp_path), "kinetics", "families")
        os.makedirs(families_root)
        _clear_gate_caches()
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        monkeypatch.setattr(rmgpy.data.rmg, "database", None, raising=False)

        self._unlistable(monkeypatch, families_root)
        try:
            with pytest.raises(getattr(module, "QuarantineEnumerationError", DatabaseError)):
                list(module.iter_quarantines())
        finally:
            _clear_gate_caches()

    def test_an_unanswered_loaded_family_without_a_quarantine_is_incomplete(
            self, tmp_path, monkeypatch):
        """
        Round 113, the reviewer's reproduction. A loaded family whose cached `quarantine` is
        None, and whose label names no directory here: `resolve_quarantine` answers
        ``(None, False)``. At `fc60e5ba4` that was dropped because the attribute was None,
        and the disk pass skipped the label as already seen -- ``([], [])``.
        """
        from rmgpy.data.kinetics import quarantine as module

        families_root = os.path.join(str(tmp_path), "kinetics", "families")
        os.makedirs(families_root)
        _clear_gate_caches()
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        _register_families(monkeypatch, {"F": None})
        monkeypatch.setattr(module, "_family_directory", lambda root, label: None)
        real_listdir = os.listdir
        monkeypatch.setattr(os, "listdir", lambda path=".": ["F"] if os.fspath(path) == families_root
                            else real_listdir(path))

        assert resolve_quarantine("F") == (None, False)
        found, incomplete = module._enumerate_quarantines()
        assert found == [] and incomplete, (
            "an unanswered family was enumerated as clean: ({0}, {1})".format(found, incomplete))
        with pytest.raises(module.QuarantineEnumerationError):
            list(module.iter_quarantines())
        _clear_gate_caches()

    def test_a_loaded_family_without_a_quarantine_behind_a_denied_directory_is_incomplete(
            self, tmp_path, monkeypatch):
        """
        The same hole without patching the module: the loaded family's own directory
        cannot be examined, so whether a manifest sits there is unknown.
        """
        from rmgpy.data.kinetics import quarantine as module

        family_path = os.path.join(str(tmp_path), "kinetics", "families", "F")
        os.makedirs(family_path)
        write_manifest(family_path)
        _clear_gate_caches()
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        _register_families(monkeypatch, {"F": None})
        _deny_access(monkeypatch, family_path)

        with pytest.raises(module.QuarantineEnumerationError):
            list(module.iter_quarantines())
        _clear_gate_caches()

    def test_an_unlistable_families_directory_keeps_the_unanswered_warning_on(
            self, tmp_path, monkeypatch):
        """
        `_database_has_any_quarantine` is the suppression test for `_warn_unanswered`.
        Behavioural at `71ae97bd5`: False, so the warning was switched off by the failure
        of the check that decides whether it is needed.
        """
        import rmgpy.data.rmg
        from rmgpy.data.kinetics.quarantine import _database_has_any_quarantine

        families_root = os.path.join(str(tmp_path), "kinetics", "families")
        os.makedirs(os.path.join(families_root, "A_Family_On_Disk"))
        write_manifest(os.path.join(families_root, "A_Family_On_Disk"))
        _clear_gate_caches()
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        monkeypatch.setattr(rmgpy.data.rmg, "database", None, raising=False)

        self._unlistable(monkeypatch, families_root)
        try:
            answer = _database_has_any_quarantine()
        finally:
            _clear_gate_caches()

        assert answer is True, "could not look was answered as nothing is quarantined"

    def test_controls_a_readable_empty_database_is_still_silent(self, tmp_path, monkeypatch):
        """The promise the module rests on: no manifest anywhere, nothing said, no raise."""
        import rmgpy.data.rmg
        from rmgpy.data.kinetics.quarantine import (_database_has_any_quarantine,
                                                    iter_quarantines)

        os.makedirs(os.path.join(str(tmp_path), "kinetics", "families", "Ordinary"))
        _clear_gate_caches()
        monkeypatch.setitem(settings, "database.directory", str(tmp_path))
        monkeypatch.setattr(rmgpy.data.rmg, "database", None, raising=False)
        messages = self._capture_warnings(monkeypatch)

        check_quarantine(self._unattributed(make_marcus()), "test")
        found = list(iter_quarantines())
        any_quarantine = _database_has_any_quarantine()
        _clear_gate_caches()

        assert messages == [] and found == [] and any_quarantine is False


class TestTheSubclassesCopyAsThemselves:
    """
    Round 111's second MEDIUM. Round 110 gave `TemplateReaction` and `LibraryReaction` a
    `copy()` that reproduces the whole object and stopped there; `DepositoryReaction` and
    `PDepReaction` kept the inherited method, which builds a base `Reaction`.
    `rmgpy/tools/isotopes.py:394` copies whatever `Reaction` it is handed.
    """

    @staticmethod
    def _depository_reaction():
        from rmgpy.data.kinetics.depository import DepositoryReaction, KineticsDepository

        reaction = DepositoryReaction(
            reactants=[Species(label="Lip", molecule=[Molecule(smiles="[Li+]")])],
            products=[Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")])],
            depository=KineticsDepository(label="Some_Family/training"),
            family="Some_Family",
            entry=Entry(index=7, label="a training reaction"))
        reaction.allow_max_rate_violation = True
        reaction.rank = 3
        reaction.comment = "a comment"
        reaction.label = "a label"
        return reaction

    @staticmethod
    def _pdep_reaction():
        from rmgpy.rmg.pdep import PDepReaction

        reaction = PDepReaction(
            reactants=[Species(label="Lip", molecule=[Molecule(smiles="[Li+]")])],
            products=[Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")])],
            network="a network stand-in")
        reaction.allow_max_rate_violation = True
        reaction.rank = 3
        reaction.comment = "a comment"
        reaction.elementary_high_p = True
        return reaction

    def test_a_depository_reaction_copies_as_one(self):
        """
        Behavioural at `13e3227b2`: the copy was a base `Reaction` and `entry` was gone --
        and `entry` is the carrier the quarantine gate reads authorship from, so a copied
        training reaction became unattributable.
        """
        from rmgpy.data.kinetics.depository import DepositoryReaction

        reaction = self._depository_reaction()
        copied = reaction.copy()

        assert isinstance(copied, DepositoryReaction)
        assert copied.entry is reaction.entry
        assert copied.depository is reaction.depository
        assert copied.family == reaction.family
        assert copied.allow_max_rate_violation is True
        assert copied.rank == 3

    def test_a_pdep_reaction_copies_as_one(self):
        """
        Behavioural at `13e3227b2`: the copy was a base `Reaction` and `network` was gone,
        which is the only thing `get_source` has to report.
        """
        from rmgpy.rmg.pdep import PDepReaction

        reaction = self._pdep_reaction()
        copied = reaction.copy()

        assert isinstance(copied, PDepReaction)
        assert copied.network is reaction.network
        assert copied.get_source() == reaction.get_source()
        assert copied.allow_max_rate_violation is True
        assert copied.elementary_high_p is True

    def test_both_subclasses_survive_a_pickle_round_trip_whole(self):
        """
        Behavioural at `13e3227b2`: the hand-written reducers turned
        ``allow_max_rate_violation=True`` into ``False`` and ``rank=3`` into ``None``.
        Round 105's fields, in a third and fourth hand-written list.
        """
        for reaction in (self._depository_reaction(), self._pdep_reaction()):
            after = pickle.loads(pickle.dumps(reaction))
            assert type(after) is type(reaction)
            assert after.allow_max_rate_violation is True
            assert after.rank == 3
            assert after.comment == "a comment"

    def test_both_subclasses_carry_their_nested_state_through_every_rmg_transport(self):
        """
        The test above checks scalars only, so nested loss could not fail it (round 112
        review). Here the payload is inside the reactants: atom ``id`` and ``props``, which
        drive resonance correspondence and group matching. Measured at `c44390545`: all
        three RMG transports carry them. A plain `pickle.dumps` does not -- every id ``-1``,
        every ``props`` empty -- and that residue is `Atom`'s own reducer in
        `rmgpy/molecule/`, outside this round's gates; it is asserted last so that an
        upstream fix announces itself here.
        """
        from multiprocessing.reduction import ForkingPickler
        from rmgpy.data.kinetics.family import complete_round_trip

        transports = {
            "copy()": lambda reaction: reaction.copy(),
            "complete_round_trip": complete_round_trip,
            "ForkingPickler": lambda reaction: pickle.loads(ForkingPickler.dumps(reaction)),
        }
        for reaction in (self._depository_reaction(), self._pdep_reaction()):
            for side in (reaction.reactants, reaction.products):
                for n, atom in enumerate(side[0].molecule[0].atoms):
                    atom.id = 100 + n
                    atom.props = {"inRing": bool(n % 2), "marker": n}
            expected = [[(atom.id, dict(atom.props)) for atom in side[0].molecule[0].atoms]
                        for side in (reaction.reactants, reaction.products)]

            for name, transport in transports.items():
                after = transport(reaction)
                assert type(after) is type(reaction), name
                got = [[(atom.id, dict(atom.props)) for atom in side[0].molecule[0].atoms]
                       for side in (after.reactants, after.products)]
                assert got == expected, name

            plain = pickle.loads(pickle.dumps(reaction))
            assert {atom.id for atom in plain.products[0].molecule[0].atoms} == {-1}, (
                "plain pickle now carries atom ids: rmgpy/molecule/ changed upstream, so "
                "retire the residue note in family.complete_round_trip")

    def test_a_deep_copy_of_either_keeps_its_class(self):
        """`deepcopy` and `copy()` are one mechanism, for these two as for the other two."""
        from copy import deepcopy

        for reaction in (self._depository_reaction(), self._pdep_reaction()):
            assert type(deepcopy(reaction)) is type(reaction)

    def test_the_two_new_fields_are_classified_explicitly(self):
        """
        The partition has no silent default, so a subclass field must be named before its
        reaction can be copied at all. Both are carried by reference, and both say why.

        **Structural** at `13e3227b2`.
        """
        from rmgpy.data.kinetics.family import _COPIED_BY_REFERENCE

        for field in ("depository", "network"):
            assert field in _COPIED_BY_REFERENCE
            assert len(_COPIED_BY_REFERENCE[field]) > 40, (
                "{0!r} is classified without a reason".format(field))


class TestWhatAnUnlistedLossyClassCosts:
    """
    The census's own question: the set of lossy classes is hand-enumerated, so what happens
    to a class nobody has added to it?

    Round 111 answered by measurement -- **quietly wrong** -- and pinned that answer in a
    test, which recorded the defect rather than closing it. Round 112 decides it:

    * **A subclass of a registered class is refused, loudly, at pickling time.** The
      registration is exact-type, so a subclass would otherwise fall through to the lossy
      reducer it inherits -- the one that names its *parent* as the class to rebuild. This
      is not hypothetical: `arkane.encorr.data.Molecule` is such a subclass, and a plain
      pickle of it comes back a base `Molecule` without its ``id``.
    * **An unrelated class with its own lossy ``__reduce__`` stays undetectable**, and that
      is accepted, with the reason: a reducer that omits a field is indistinguishable from
      a class that does not have one. What keeps the table honest there is the census,
      which round-trips every class reachable in a reaction's state.
    * **Private Cython memos are accepted as caches** -- see
      `test_the_private_cython_state_is_exactly_the_three_memos`.
    """

    @staticmethod
    def _subclass_of(base):
        """An unregistered subclass of `base` that carries one field of its own."""
        return type("Unregistered" + base.__name__, (base,), {})

    def test_an_unregistered_subclass_is_refused_by_the_copy_transport(self):
        """
        Behavioural at `71ae97bd5`: `complete_round_trip` returned a base `Molecule`,
        class and extra field gone, without a word.
        """
        from rmgpy.data.kinetics.family import complete_round_trip

        Sub = self._subclass_of(Molecule)
        molecule = Sub(smiles="CC")
        molecule.extra = "carried"

        with pytest.raises(pickle.PicklingError, match="UnregisteredMolecule"):
            complete_round_trip(molecule)
        # Nested is the case that matters: the subclass sits inside a reaction's state.
        with pytest.raises(pickle.PicklingError, match="UnregisteredMolecule"):
            complete_round_trip([Species(molecule=[molecule])])

    def test_an_unregistered_atom_subclass_is_refused_too(self):
        """`Atom` is registered with the *other* reducer; the refusal must not care which."""
        from rmgpy.data.kinetics.family import complete_round_trip
        from rmgpy.molecule.molecule import Atom

        with pytest.raises(pickle.PicklingError, match="UnregisteredAtom"):
            complete_round_trip(self._subclass_of(Atom)(element="C"))

    def test_the_multiprocessing_pickler_refuses_a_subclass_loaded_before_install(self):
        """
        Behavioural at `71ae97bd5`: ``ForkingPickler`` fell through to the inherited
        reducer. Round 112 registered a refusal per subclass loaded when
        `install_complete_reducers` ran; round 113 resolves by MRO at reduce time instead,
        which covers this case and the one defined after it (next test).
        """
        from multiprocessing.reduction import ForkingPickler
        from rmgpy.data.kinetics.family import install_complete_reducers

        Sub = self._subclass_of(Molecule)
        install_complete_reducers()
        try:
            with pytest.raises(pickle.PicklingError, match="UnregisteredMolecule"):
                ForkingPickler.dumps([Sub(smiles="C")])
        finally:
            ForkingPickler._extra_reducers.pop(Sub, None)

    def test_the_multiprocessing_pickler_refuses_a_subclass_defined_after_install(self):
        """
        Round 113. `install_complete_reducers` ran at import, long before this class
        existed, and nothing re-runs it. Behavioural at `fc60e5ba4`: ``ForkingPickler``
        rebuilt a base `Molecule` and the extra field was gone.
        """
        from multiprocessing.connection import Pipe
        from multiprocessing.reduction import ForkingPickler

        Sub = self._subclass_of(Molecule)
        assert Sub not in ForkingPickler._extra_reducers
        molecule = Sub(smiles="C")
        molecule.extra = "carried"

        with pytest.raises(pickle.PicklingError, match="UnregisteredMolecule"):
            ForkingPickler.dumps([Species(molecule=[molecule])])
        # The transport itself, not only its pickler: `Connection.send` is what `Pool` uses.
        sender, receiver = Pipe()
        try:
            with pytest.raises(pickle.PicklingError, match="UnregisteredMolecule"):
                sender.send(molecule)
        finally:
            sender.close()
            receiver.close()

    def test_a_verdict_is_not_inherited_through_an_equal_comparing_metaclass(self):
        """
        Round 113 rework. The verdict cache was a `WeakKeyDictionary`, whose lookup goes
        through the metaclass's ``__eq__``/``__hash__`` rather than class identity. A benign
        class cached ``False``; a late lossy `Molecule` subclass whose metaclass compares
        equal to it then read that ``False`` and crossed the pickler as a base `Molecule`.
        Behavioural at `2f35a866a`.
        """
        from multiprocessing.reduction import ForkingPickler

        class Colliding(type):
            def __eq__(cls, other):
                return isinstance(other, Colliding)

            def __hash__(cls):
                return 0

        class Benign(metaclass=Colliding):
            pass

        # The verdict is cached before pickle finds a local class unpicklable by reference.
        with pytest.raises((AttributeError, pickle.PicklingError), match="local object"):
            ForkingPickler.dumps(Benign())
        Sub = Colliding("UnregisteredMolecule", (Molecule,), {})
        molecule = Sub(smiles="C")
        molecule.extra = "carried"

        with pytest.raises(pickle.PicklingError, match="UnregisteredMolecule"):
            ForkingPickler.dumps(molecule)

    def test_registered_classes_still_travel(self):
        """The control: the refusal is for subclasses, not for the classes themselves."""
        from multiprocessing.reduction import ForkingPickler
        from rmgpy.data.kinetics.family import complete_round_trip
        from rmgpy.molecule.fragment import Fragment

        for value in (Molecule(smiles="CC"), Species(molecule=[Molecule(smiles="C")]),
                      Fragment().from_smiles_like_string("CCR")):
            assert type(complete_round_trip(value)) is type(value)
            assert type(pickle.loads(ForkingPickler.dumps(value))) is type(value)

    def test_every_subclass_in_the_tree_is_registered_or_accepted(self):
        """
        Mechanical: import every module of `rmgpy` and `arkane` in a fresh interpreter,
        walk ``__subclasses__()`` recursively under each registered class, and fail on a
        subclass that is neither registered nor listed below with a reason. A subprocess,
        because importing the whole tree into the test process would leak module state
        into every later test.
        """
        from rmgpy.data.kinetics.family import _ACCEPTED_UNREGISTERED_SUBCLASSES

        script = textwrap.dedent("""
            import importlib, pkgutil, logging
            logging.disable(logging.CRITICAL)
            import rmgpy, arkane
            from rmgpy.data.kinetics.family import _LOSSY_REDUCERS
            for package in (rmgpy, arkane):
                for module in pkgutil.walk_packages(package.__path__, package.__name__ + "."):
                    try:
                        importlib.import_module(module.name)
                    except Exception:
                        pass
            def walk(cls):
                for sub in cls.__subclasses__():
                    yield sub
                    yield from walk(sub)
            for base in _LOSSY_REDUCERS:
                for sub in walk(base):
                    if sub not in _LOSSY_REDUCERS:
                        print(sub.__module__ + "." + sub.__qualname__)
        """)
        result = subprocess.run([sys.executable, "-c", script], capture_output=True,
                                text=True, timeout=600)
        assert result.returncode == 0, result.stderr
        found = set(result.stdout.split())
        assert found == set(_ACCEPTED_UNREGISTERED_SUBCLASSES), (
            "unregistered subclasses of a lossy class: {0}".format(sorted(found)))
        for name, reason in _ACCEPTED_UNREGISTERED_SUBCLASSES.items():
            assert len(reason) > 40, "{0} is accepted without a reason".format(name)

    def test_the_private_cython_state_is_exactly_the_three_memos(self):
        """
        `writable_fields` cannot see a ``cdef`` attribute that is neither ``public`` nor
        ``readonly``, so no transport carries one. Accepted, for a reason pinned here: the
        only such state on the registered classes is ``_fingerprint``, ``_inchi`` and
        ``_smiles`` on `Molecule` and `Species`, each a memo of a value derived from the
        graph and recomputed when absent. Parsed from the ``.pxd`` so that a new private
        field fails this test instead of vanishing from every copy.
        """
        import rmgpy

        root = os.path.dirname(rmgpy.__file__)
        private = {}
        for pxd in ("molecule/graph.pxd", "molecule/molecule.pxd", "molecule/fragment.pxd",
                    "species.pxd"):
            path = os.path.join(root, pxd)
            if not os.path.exists(path):
                continue
            current = None
            for line in open(path):
                klass = re.match(r"cdef class (\w+)", line)
                if klass:
                    current = klass.group(1)
                    continue
                field = re.match(r"    cdef (?!public\b|readonly\b)\w[\w\.\[\]]*\s+(\w+)\s*$",
                                 line)
                if field and current:
                    private.setdefault(current, set()).add(field.group(1))
        memos = {"_fingerprint", "_inchi", "_smiles"}
        assert private == {"Molecule": memos, "Species": memos}

    def test_the_memos_are_recomputed_after_a_copy(self):
        """The half of the acceptance that is behaviour: a copy answers as the original."""
        from rmgpy.data.kinetics.family import complete_round_trip

        molecule = Molecule(smiles="C=CC")
        species = Species(molecule=[Molecule(smiles="C=CC")])
        expected = (molecule.fingerprint, molecule.smiles, molecule.to_inchi())
        expected_species = (species.fingerprint, species.smiles, species.inchi)

        copied, copied_species = complete_round_trip([molecule, species])
        assert (copied.fingerprint, copied.smiles, copied.to_inchi()) == expected
        assert (copied_species.fingerprint, copied_species.smiles,
                copied_species.inchi) == expected_species

    def test_the_only_measured_member_of_the_table_is_atom(self):
        """
        What "hand-enumerated" means precisely. `Atom`'s losses are *measured* by
        `fields_the_reducer_drops`, on a probe, so upstream can add or remove a field and
        the restore follows. The other four are judgements written down with reasons, and
        the census test is what keeps them honest.

        **Structural** at `13e3227b2` in part -- `COMPLETE_REDUCERS` is new -- and a
        statement of scope rather than a defect.
        """
        from rmgpy.data.kinetics.family import _LOSSY_REDUCERS, _REDUCER_PROBES
        from rmgpy.molecule.molecule import Atom

        assert set(_REDUCER_PROBES) == {Atom}
        assert set(_REDUCER_PROBES) <= set(_LOSSY_REDUCERS)
        for _cls, (_reducer, reason) in _LOSSY_REDUCERS.items():
            assert len(reason) > 40, "a table entry carries no reason"
