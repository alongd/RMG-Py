"""I-059 -- the producer writes EVIDENCE for ``chain_mass_defect_g_mol``.

Every spawn path that creates or grows a per-chain mass defect performs the
same addition, ``defect_dst = defect_src + shed``, and until I-059 it threw
both operands away: the sidecar carried the RESULT and the spawn CHANNEL, but
nothing a different program in a different environment could use to check the
NUMBER. A channel name cannot tell 1.008 g/mol from 16, so the downstream
consumer could only census-warn and apply the defect-aware condensed-mass
formula on trust.

The producer now records the addition itself, next to its result::

    "chain_mass_defect_provenance": {
        "source": ..., "parent_pool": ...,
        "inherited_g_mol": ..., "shed_g_mol": ...
    }

with ``inherited_g_mol + shed_g_mol == chain_mass_defect_g_mol`` exactly.
This module pins that the block is emitted where a defect is, that it travels
with the defect through ``Polymer.copy()``, that the serializer refuses to
write a block that has drifted from the value it claims to explain, and --
the point of the ticket -- that the defect VALUE is untouched.
"""

import pytest

from rmgpy.molecule import Molecule
from rmgpy.polymer import Polymer, _serialize_pool_for_sidecar

PROV = "chain_mass_defect_provenance"
DEFECT = "chain_mass_defect_g_mol"

#: The real novolak repeat unit and one of its H-loss features, lifted from a
#: completed run's sidecar (poly_105, rmg_commit b34c787ec) so the pinned
#: value below is the one the field actually produces, not a fixture's
#: invention.
NOVOLAK_MONOMER = "[CH2]C1=CC=C(C)C([CH2])=C1O"
NOVOLAK_FEATURE = "[CH2]C1=CC=C(C)C([CH2])=C1[O]"
#: RMG's mass for one hydrogen atom -- the exact float that run serialized.
H_G_MOL = 1.0079710045829415


@pytest.fixture
def novolak():
    return Polymer(
        label="novolak",
        monomer=NOVOLAK_MONOMER,
        end_groups=["[H]", "[H]"],
        cutoff=3,
        Mn=1000.0,
        Mw=3000.0,
        initial_mass=99.0,
    )


@pytest.fixture
def h_loss_daughter(novolak):
    feat = Molecule().from_smiles(NOVOLAK_FEATURE)
    return novolak._born_at_zero_mod_daughter(
        feat, source="radical_feature_h_loss")


class TestDefectProvenanceIsWritten:
    def test_an_ordinary_pool_has_neither_defect_nor_evidence(self, novolak):
        assert novolak.chain_mass_defect_g_mol == 0.0
        assert novolak.chain_mass_defect_provenance is None
        assert PROV not in _serialize_pool_for_sidecar(novolak)

    def test_h_loss_daughter_carries_the_evidence(self, h_loss_daughter):
        assert h_loss_daughter.chain_mass_defect_provenance == {
            "source": "radical_feature_h_loss",
            "parent_pool": "novolak",
            "inherited_g_mol": 0.0,
            "shed_g_mol": H_G_MOL,
        }

    def test_the_defect_value_is_unchanged(self, h_loss_daughter):
        """The ticket is provenance, not arithmetic: the number must be the
        one a pre-I-059 run wrote for this exact pool."""
        assert h_loss_daughter.chain_mass_defect_g_mol == H_G_MOL

    def test_the_block_reaches_the_sidecar_beside_the_value(
            self, h_loss_daughter):
        d = _serialize_pool_for_sidecar(h_loss_daughter)
        assert d[DEFECT] == H_G_MOL
        assert d[PROV]["shed_g_mol"] == H_G_MOL
        assert d[PROV]["parent_pool"] == "novolak"

    def test_the_arithmetic_closes_exactly(self, h_loss_daughter):
        """Not "to within a tolerance" -- the consumer re-adds the same two
        IEEE doubles the producer added, so the sum is bit-for-bit equal."""
        d = _serialize_pool_for_sidecar(h_loss_daughter)
        assert (d[PROV]["inherited_g_mol"] + d[PROV]["shed_g_mol"]
                == d[DEFECT])

    def test_a_second_generation_daughter_inherits_and_sheds(
            self, h_loss_daughter):
        feat = Molecule().from_smiles(NOVOLAK_FEATURE)
        gen2 = h_loss_daughter._born_at_zero_mod_daughter(
            feat, source="radical_feature_h_loss")
        prov = gen2.chain_mass_defect_provenance
        assert prov["parent_pool"] == h_loss_daughter.label
        assert prov["inherited_g_mol"] == H_G_MOL
        assert (prov["inherited_g_mol"] + prov["shed_g_mol"]
                == gen2.chain_mass_defect_g_mol)


class TestEvidenceTravelsWithTheValue:
    def test_copy_carries_the_block(self, h_loss_daughter):
        """A copy that kept the number and dropped its provenance is exactly
        the unverifiable shape this ticket exists to remove."""
        other = h_loss_daughter.copy(deep=True)
        assert other.chain_mass_defect_g_mol == H_G_MOL
        assert other.chain_mass_defect_provenance == \
            h_loss_daughter.chain_mass_defect_provenance

    def test_copy_does_not_alias_the_block(self, h_loss_daughter):
        other = h_loss_daughter.copy(deep=True)
        other.chain_mass_defect_provenance["shed_g_mol"] = 99.0
        assert h_loss_daughter.chain_mass_defect_provenance["shed_g_mol"] \
            == H_G_MOL

    def test_a_serialized_copy_is_accepted_by_the_loader(
            self, novolak, h_loss_daughter):
        """The copy shape is the trap here. ``Polymer.copy()`` deliberately
        does NOT carry ``spawn_metadata``, so a copied defect-bearing pool
        serializes with the ``{'source': 'input'}`` absence sentinel next to
        a block naming ``radical_feature_h_loss``. A loader that demanded
        those two agree would hard-reject a real producer artifact."""
        from rmgpy.tools.polymer_moments_runner import (
            validate_chain_mass_defect_provenance)
        copied = h_loss_daughter.copy(deep=True)
        entry = _serialize_pool_for_sidecar(copied)
        assert entry["spawn_event_metadata"] == {"source": "input"}
        assert entry[PROV]["source"] == "radical_feature_h_loss"
        parent = _serialize_pool_for_sidecar(novolak)
        confirmed, unresolved = validate_chain_mass_defect_provenance(
            [parent, entry])
        assert confirmed == {entry["label"]}
        assert unresolved == {}

    def test_a_copy_naming_a_different_channel_is_still_refused(
            self, novolak, h_loss_daughter):
        """The absence sentinel widens the cross-pin, it does not remove it:
        a pool that DOES record its own spawn event must still agree."""
        from rmgpy.tools.polymer_moments_runner import (
            validate_chain_mass_defect_provenance)
        entry = _serialize_pool_for_sidecar(h_loss_daughter)
        entry[PROV]["source"] = "concerted_loss_ejection"
        with pytest.raises(ValueError, match="contradicts"):
            validate_chain_mass_defect_provenance(
                [_serialize_pool_for_sidecar(novolak), entry])


class TestSerializerRefusesDriftedEvidence:
    def test_drifted_block_is_refused_at_write_time(self, h_loss_daughter):
        """The write is the reconciliation point: a defect mutated after
        spawn without updating its evidence must not be laundered into a
        consumer's mass balance as unverifiable provenance."""
        h_loss_daughter.chain_mass_defect_g_mol = 2.0 * H_G_MOL
        with pytest.raises(ValueError, match="does not reconstruct"):
            _serialize_pool_for_sidecar(h_loss_daughter)

    def test_a_defect_with_no_block_still_serializes(self, h_loss_daughter):
        """Absence stays legal -- pools whose defect predates the evidence
        path must keep serializing, or every old code path breaks."""
        h_loss_daughter.chain_mass_defect_provenance = None
        d = _serialize_pool_for_sidecar(h_loss_daughter)
        assert d[DEFECT] == H_G_MOL
        assert PROV not in d


class TestReferenceLoaderMirror:
    """RMG's own reference loader validates the evidence with the same rules
    TA does (the producer/consumer mirror property). If these two ever
    disagree, one of them is wrong about what the artifact means."""

    @staticmethod
    def _pools(**overrides):
        parent = {"label": "novolak",
                  "spawn_event_metadata": {"source": "input"}}
        child = {"label": "novolak_mod",
                 "spawn_event_metadata": {"source": "radical_feature_h_loss"},
                 DEFECT: H_G_MOL,
                 PROV: {"source": "radical_feature_h_loss",
                        "parent_pool": "novolak",
                        "inherited_g_mol": 0.0,
                        "shed_g_mol": H_G_MOL}}
        child[PROV].update(overrides)
        return [parent, child]

    def test_well_formed_evidence_is_confirmed(self):
        from rmgpy.tools.polymer_moments_runner import (
            validate_chain_mass_defect_provenance)
        confirmed, unresolved = validate_chain_mass_defect_provenance(
            self._pools())
        assert confirmed == {"novolak_mod"}
        assert unresolved == {}

    def test_arithmetic_that_does_not_close_is_refused(self):
        from rmgpy.tools.polymer_moments_runner import (
            validate_chain_mass_defect_provenance)
        with pytest.raises(ValueError, match="does not reconstruct"):
            validate_chain_mass_defect_provenance(
                self._pools(shed_g_mol=2.0 * H_G_MOL))

    def test_source_must_match_the_spawn_event(self):
        from rmgpy.tools.polymer_moments_runner import (
            validate_chain_mass_defect_provenance)
        with pytest.raises(ValueError, match="contradicts"):
            validate_chain_mass_defect_provenance(
                self._pools(source="side_group_homolysis"))

    def test_missing_parent_is_reported_not_confirmed(self):
        from rmgpy.tools.polymer_moments_runner import (
            validate_chain_mass_defect_provenance)
        confirmed, unresolved = validate_chain_mass_defect_provenance(
            self._pools(parent_pool="gone"))
        assert confirmed == set()
        assert unresolved == {"novolak_mod": "gone"}

    def test_absence_is_legal(self):
        from rmgpy.tools.polymer_moments_runner import (
            validate_chain_mass_defect_provenance)
        pools = self._pools()
        del pools[1][PROV]
        assert validate_chain_mass_defect_provenance(pools) == (set(), {})
