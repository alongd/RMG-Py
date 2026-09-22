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
"""

import inspect
import logging
import os

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
)
from rmgpy.exceptions import DatabaseError, QuarantinedKineticsError
from rmgpy.kinetics.arrhenius import Arrhenius, Marcus
from rmgpy.molecule import Molecule
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species

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
def registered(monkeypatch, quarantine):
    """
    Register the synthetic quarantined family in the kinetics database singleton.

    `check_quarantine` resolves a reaction's family label the same way the rest of RMG
    does, so this is what makes the gate see the synthetic family.
    """

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

    for name in ("_UNATTRIBUTED_WARNED", "_UNANSWERED_WARNED", "_DISK_QUARANTINE_CACHE"):
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
        assert all(isinstance(e.data, Marcus) for e in affected["rules"] + affected["training"])

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
        assert entry.data.lmbd_i_coefs.value_si[0] == pytest.approx(51487.7)
        assert entry.data.lmbd_o.value_si == pytest.approx(0.0)


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
