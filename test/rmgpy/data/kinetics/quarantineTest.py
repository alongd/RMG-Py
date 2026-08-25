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
import os

import pytest

from rmgpy import settings
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.quarantine import (
    QUARANTINE_FILENAME,
    KineticsQuarantine,
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
