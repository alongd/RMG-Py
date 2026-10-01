#!/usr/bin/env python
"""I-067 round 2: invariants of the measurement instruments themselves.

These cover the round-2 review findings that do not need a full RMG generation to
demonstrate. The ones that DO need a run (that the time buckets sum to the enlargement
wall clock) are demonstrated by the residual line the reporter now emits per enlargement.

Every test here was written to FAIL against the round-1 code and was confirmed red before
the fix landed; see the contract's Evidence section.
"""
import logging

import pytest

from rmgpy.constraints import _record_refusal, get_generation_census, reset_generation_census
from rmgpy.data.kinetics import family as family_mod
from rmgpy.molecule import Molecule
from rmgpy.rmg.model import CoreEdgeReactionModel


class TestWastedBuildProfileReset:
    """A module-global accumulator must not leak across runs in one process."""

    def test_reset_zeroes_every_field(self):
        p = family_mod.get_wasted_build_profile()
        for key in p:
            p[key] = 7 if isinstance(p[key], int) else 7.0

        family_mod.reset_wasted_build_profile()

        p = family_mod.get_wasted_build_profile()
        assert p["procnum"] == 1
        for key, value in p.items():
            if key == "procnum":
                continue
            assert value == 0, "{0} survived the reset with {1!r}".format(key, value)

    def test_reset_keeps_the_same_dict_object(self):
        """Callers hold the dict returned by get_wasted_build_profile(); rebinding the
        module global instead of clearing in place would leave them on a stale copy."""
        before = family_mod.get_wasted_build_profile()
        family_mod.reset_wasted_build_profile()
        assert family_mod.get_wasted_build_profile() is before


class TestPartitionIsComplete:
    """The *_total_s buckets must cover every exit from the timed region."""

    def test_reporter_and_accumulator_agree_on_the_bucket_set(self):
        acc = {k for k in family_mod.get_wasted_build_profile() if k.endswith("_total_s")}
        assert acc == set(CoreEdgeReactionModel.WASTED_BUILD_BUCKETS), (
            "the reporter sums a different set of buckets than the accumulator keeps, so "
            "its 'attributed' total would silently omit one"
        )

    def test_every_charge_site_names_a_real_bucket(self):
        """A typo'd bucket name would create a new key rather than raising, and the
        reporter would then never see that time."""
        import inspect
        import re

        src = inspect.getsource(family_mod.KineticsFamily._generate_product_structures)
        charged = set(re.findall(r"_charge\(['\"]([a-z_]+)['\"]\)", src))
        assert charged, "found no _charge() call sites -- the extractor is broken"
        unknown = charged - set(CoreEdgeReactionModel.WASTED_BUILD_BUCKETS)
        assert not unknown, "these _charge() sites name buckets nobody sums: {0}".format(unknown)


class TestProcnumGuard:
    """Under forked generation the accumulator is meaningless and must say so."""

    def _emit(self, caplog, procnum):
        family_mod.reset_wasted_build_profile()
        p = family_mod.get_wasted_build_profile()
        p["procnum"] = procnum
        p["refused_recipe_s"] = 5.0
        p["refused_total_s"] = 5.0
        p["accepted_total_s"] = 1.0
        model = CoreEdgeReactionModel.__new__(CoreEdgeReactionModel)
        with caplog.at_level(logging.INFO):
            model.log_wasted_build_profile(wall_seconds=10.0)
        return caplog.text

    def test_multiprocess_profile_reports_no_fraction(self, caplog):
        text = self._emit(caplog, procnum=8)
        assert "NOT VALID" in text
        assert "procnum=8" in text
        assert "%" not in text, "a share of wall time was printed for a forked run"

    def test_single_process_profile_reports_the_fraction(self, caplog):
        text = self._emit(caplog, procnum=1)
        assert "NOT VALID" not in text
        assert "50.0%" in text, "the refused share of a 5 s / 10 s enlargement should be 50%"

    def test_note_generation_procnum_keeps_the_maximum(self):
        family_mod.reset_wasted_build_profile()
        original = family_mod.PROFILE_WASTED_BUILDS
        family_mod.PROFILE_WASTED_BUILDS = True
        try:
            family_mod.note_generation_procnum(4)
            family_mod.note_generation_procnum(1)
        finally:
            family_mod.PROFILE_WASTED_BUILDS = original
        assert family_mod.get_wasted_build_profile()["procnum"] == 4, (
            "one parallel enlargement taints the whole cumulative profile, so the guard "
            "must latch rather than track the most recent call"
        )


class TestCensusDoesNotMutateWhatItMeasures:
    """to_smiles() reorders mol.atoms in place; the census must render a copy."""

    @staticmethod
    def _snapshot(mol):
        return [
            (a.element.symbol, a.radical_electrons, a.lone_pairs, a.charge,
             getattr(a, "sorting_label", None))
            for a in mol.atoms
        ]

    @pytest.mark.parametrize("smiles", ["Oc1ccccc1", "Cc1ccccc1O", "[CH2]CCC[CH2]", "[CH]C"])
    def test_record_refusal_leaves_the_structure_untouched(self, smiles):
        reset_generation_census()
        mol = Molecule(smiles=smiles)
        before = self._snapshot(mol)

        _record_refusal(mol, "gas", "Exceeded maximumCarbonAtoms: 30")

        assert self._snapshot(mol) == before, (
            "the refusal census permuted the structure it was measuring"
        )
        assert sum(get_generation_census()["refused"].values()) == 1

    def test_the_refusal_is_still_recorded_with_its_smiles(self):
        reset_generation_census()
        mol = Molecule(smiles="Oc1ccccc1")
        _record_refusal(mol, "gas", "Exceeded maximumCarbonAtoms: 30")
        bucket = get_generation_census()["smiles"]["gas"]
        assert len(bucket) == 1
        (count, heavy, carbon), = bucket.values()
        assert (count, heavy, carbon) == (1, 7, 6)


class TestCapAccounting:
    """Past the cap the census drops repeats too, and must not hide that."""

    def test_dropped_refusals_are_counted(self, monkeypatch):
        monkeypatch.setattr("rmgpy.constraints._CENSUS_SMILES_CAP", 1)
        reset_generation_census()
        _record_refusal(Molecule(smiles="Oc1ccccc1"), "gas", "Exceeded maximumCarbonAtoms: 30")
        for _ in range(3):
            _record_refusal(Molecule(smiles="Cc1ccccc1O"), "gas",
                            "Exceeded maximumCarbonAtoms: 30")

        census = get_generation_census()
        assert census["smiles_capped"] is True
        assert census["smiles_dropped"] == 3, (
            "refusals dropped after the cap must be counted, otherwise the per-structure "
            "counts read as exact when they are lower bounds"
        )
        # The exact counters are unaffected by the cap.
        assert sum(census["refused"].values()) == 4
        assert sum(census["by_heavy"]["gas"].values()) == 4
