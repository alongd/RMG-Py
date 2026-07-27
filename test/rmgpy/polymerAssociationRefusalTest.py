#!/usr/bin/env python3

"""
Pins for the class-level refusal of the ASSOCIATION orientation -- a row whose
only Polymer participant is a PRODUCT (rmgpy/polymer.py,
:func:`stamp_gas_association_refusal`).

The rule is structural, not a taste in shapes: ``stamp_polymer_flux_archetype``
is reachable only under ``if polymer_reactants:`` (rmgpy/rmg/model.py:684 and
:760), and ``restamp_flipped_polymer_archetype`` only re-derives an archetype
for a row that already carries one. So a polymer-as-product-only row can NEVER
be archetype-stamped: its only alternative to refusal is arriving unstamped at
the solver rebuild and killing the run at the r71 guard.

Grounded in two runs on 2026-07-27 -- an EPDM copolymer deck and a plain PE
homopolymer control -- which died identically on the row pinned first below.
"""

from rmgpy.molecule import Molecule
from rmgpy.polymer import Polymer, stamp_gas_association_refusal
from rmgpy.reaction import Reaction
from rmgpy.species import Species


def _pe_pool():
    """A polyolefin pool whose baseline proxy is the small saturated C7H16."""
    return Polymer(label="PE", monomer="[CH2][CH2]",
                   end_groups=["[CH3]", "[H]"], cutoff=3,
                   Mn=5000.0, Mw=10000.0, initial_mass=1.0)


def _spc(smiles):
    return Species(molecule=[Molecule().from_smiles(smiles)])


def _rxn(reactants, products):
    return Reaction(reactants=reactants, products=products, reversible=True)


class TestAssociationOrientationRefusal:

    def test_carbene_plus_closed_shell_association_is_refused(self):
        """
        THE regression case: ``[CH2] + CCCCCC <=> <pool>``. CH2 inserting into
        hexane reconstitutes the C7H16 baseline proxy of any polyolefin pool at
        cutoff=3.

        It escaped every enumerated shape: the r63 all-radicals conjunct fails
        (hexane is closed-shell) and the r82 impostor conjunct correctly
        declines (hexane is 86 g/mol against the 300 g/mol absolute chain-scale
        floor). Left unrefused, it reached the solver rebuild unstamped and
        killed the run.
        """
        rxn = _rxn([_spc("[CH2]"), _spc("CCCCCC")], [_pe_pool()])
        stamp_gas_association_refusal(rxn)
        assert rxn.polymer_refused is True
        assert rxn.polymer_refused_accumulating is False  # conduit-deferred

    def test_association_refused_regardless_of_partner_size_or_radicals(self):
        """The orientation decides; size and radical count no longer gate it."""
        for pair in (("[CH3]", "CCCCC"),      # radical + closed-shell
                     ("CC", "CCCC"),          # both closed-shell
                     ("[CH3]", "[CH2]C")):    # both radicals (legacy r63 case)
            rxn = _rxn([_spc(pair[0]), _spc(pair[1])], [_pe_pool()])
            stamp_gas_association_refusal(rxn)
            assert rxn.polymer_refused is True, \
                'association orientation left live for %s' % (pair,)

    def test_homolysis_orientation_is_not_widened(self):
        """
        Negative control: the orientation that CAN be archetype-stamped keeps
        the narrow r63 conjunct, so volatile-producing chain chemistry stays
        live. Widening this side would zero the flux of real depropagation.
        """
        rxn = _rxn([_pe_pool()], [_spc("C=C"), _spc("C=CC")])
        stamp_gas_association_refusal(rxn)
        assert rxn.polymer_refused is False

    def test_homolysis_all_radicals_still_refused(self):
        """The legacy r63 conjunct on the homolysis side is untouched."""
        rxn = _rxn([_pe_pool()], [_spc("[CH3]"), _spc("[CH2]C")])
        stamp_gas_association_refusal(rxn)
        assert rxn.polymer_refused is True

    def test_both_sides_condensed_untouched(self):
        """
        Routing chemistry (pool on both sides) is not this rule's business.

        The daughter is marked FEATURE, matching the repo's own idiom in
        ``test_impostor_refusal_leaves_feature_pool_abstraction_alone``: an
        abstraction returning the IDENTICAL proxy on both sides is the r74
        same-proxy 'H fountain' degenerate shape, which is refused for its own
        adjudicated reason and would not isolate the orientation rule.
        """
        from rmgpy.polymer import PolymerClass
        pool = _pe_pool()
        daughter = pool.copy()
        daughter._reacted_class = PolymerClass.FEATURE
        rxn = _rxn([pool, _spc("[H]")], [daughter, _spc("[H][H]")])
        stamp_gas_association_refusal(rxn)
        assert rxn.polymer_refused is False

    def test_gas_only_rows_untouched(self):
        """No condensed participant anywhere: never this rule's business."""
        rxn = _rxn([_spc("[CH2]"), _spc("CCCCCC")], [_spc("CCCCCCC")])
        stamp_gas_association_refusal(rxn)
        assert rxn.polymer_refused is False

    def test_copolymer_pool_behaves_identically(self):
        """
        The gap is orientation-driven, not copolymer-driven: an EPDM-style
        composition pool refuses the same row the same way. This is the
        assertion form of the PE control run that cleared the copolymer
        feature of causing the failure.
        """
        epdm = Polymer(
            label="EPDM",
            monomers=[{'monomer': '[CH2][CH2]', 'fraction': 0.7101},
                      {'monomer': '[CH2][CH](C)', 'fraction': 0.2761},
                      {'monomer': 'CC=C1CC2[CH][CH]C1C2', 'fraction': 0.0138}],
            end_groups=["[CH3]", "[H]"], cutoff=3,
            Mn=5000.0, Mw=10000.0, initial_mass=1.0)
        rxn = _rxn([_spc("[CH2]"), _spc("CCCCCC")], [epdm])
        stamp_gas_association_refusal(rxn)
        assert rxn.polymer_refused is True
