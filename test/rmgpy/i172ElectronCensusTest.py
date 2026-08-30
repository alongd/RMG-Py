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
The free electron is judged by charge, not by element count.

``rmgpy.molecule.element.element_list`` begins with element ``e``, so before
I-172 ``Reaction.is_balanced`` demanded that the free-electron COUNT match across
the arrow, and returned ``False`` from that comparison before reaching the charge
comparison at the end of the method.  Every ionisation, attachment and
recombination written with explicit electrons was therefore refused, however well
its charge balanced -- and since one refusal aborts a whole kinetics library, a
plasma library carrying such entries contributed no reactions at all.

The electron is now skipped in the element census and judged solely by the charge
comparison, which already accounts for it: each ``e-`` species contributes -1
through the per-atom charge sum, and the scalar ``Reaction.electrons`` is folded
in immediately above.  It stops being judged twice -- correctly as charge, and
incorrectly as a conserved element -- not once.

**The negative controls are the point of this module.**  Widening what balances
is only safe if the charge comparison still catches what the census used to.  So
every accept below is matched by a refusal that differs from it in charge alone,
and :func:`test_electron_gain_without_charge_change_is_still_refused` is the
direct falsifier: delete the ``return reactants_net_charge == products_net_charge``
line and it goes red immediately, while the accepts stay green.

There is also a structural reason the widening is narrow, pinned in
:func:`test_charge_conserving_electron_change_implies_matching_ionisation`.  A
free electron is the only participant whose element is ``e``, and it carries
charge -1 exactly.  So if the non-electron element census matches and the total
charge is conserved, a change of `n` in the free-electron count forces a
compensating charge change of `+n` on the heavy species -- which is the
definition of ionisation by `n` electrons.  The class of newly-accepted reactions
is therefore exactly "correctly written ionisation/attachment", not "anything
with a plausible-looking electron count".

No database is needed: every participant is built from an adjacency list.
"""

import pytest

from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species

################################################################################

#: Spellings match ``test/rmgpy/i108ElectronRepresentationMatrixTest.py`` and the
#: plasma database's own dictionaries.
ADJACENCY_LISTS = {
    'e': '1 e u1 p0 c-1',
    'Ar': '1 Ar u0 p4 c0',
    'Arplus': '1 Ar u1 p3 c+1',
    'N': '1 N u3 p1 c0',
    'N2': '1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}',
    'N2plus': '1 N u0 p1 c0 {2,T}\n2 N u1 p0 c+1 {1,T}',
}


def _s(label):
    """A fresh :class:`Species`. Fresh per call, never shared between a
    reactant and a product list."""
    return Species(label=label, molecule=[Molecule().from_adjacency_list(ADJACENCY_LISTS[label])])


def _rxn(reactants, products, electrons=0):
    return Reaction(reactants=[_s(label) for label in reactants],
                    products=[_s(label) for label in products],
                    electrons=electrons, reversible=False)


def _charges(reaction):
    """``(reactant charge, product charge)``, summed the way ``is_balanced``
    sums them and with the scalar ``electrons`` folded in the same way."""
    def total(species_list):
        return sum(atom.charge for spc in species_list for atom in spc.molecule[0].atoms)

    left, right = total(reaction.reactants), total(reaction.products)
    if reaction.electrons < 0:
        left += reaction.electrons
    elif reaction.electrons > 0:
        right -= reaction.electrons
    return left, right


################################################################################
# What must start balancing
################################################################################

def test_electron_impact_ionisation_balances_by_charge():
    """``Ar + e- => Ar+ + e- + e-`` conserves charge (-1 = -1) and must balance,
    although its free-electron count changes 1 -> 2.

    This is the entry the 5 torr argon deck needs, and the shape 51 of the 52
    rejected ``PlasmaAir`` entries share.  It is the assertion that is RED before
    the change.
    """
    rxn = _rxn(['Ar', 'e'], ['Arplus', 'e', 'e'])
    assert _charges(rxn) == (-1, -1)
    assert rxn.is_balanced()


def test_half_written_ionisation_balances_by_charge():
    """``Ar => Ar+ + e-``: the electron appears on one side only, and charge is
    still conserved (0 = 0)."""
    rxn = _rxn(['Ar'], ['Arplus', 'e'])
    assert _charges(rxn) == (0, 0)
    assert rxn.is_balanced()


def test_recombination_balances_by_charge():
    """``Ar+ + e- => Ar``: the attachment direction, charge 0 = 0."""
    rxn = _rxn(['Arplus', 'e'], ['Ar'])
    assert _charges(rxn) == (0, 0)
    assert rxn.is_balanced()


def test_three_body_recombination_balances_by_charge():
    """``Ar+ + e- + e- => Ar + e-``: the electron is both a reactant and the
    third body, and the count changes 2 -> 1 with charge -1 = -1."""
    rxn = _rxn(['Arplus', 'e', 'e'], ['Ar', 'e'])
    assert _charges(rxn) == (-1, -1)
    assert rxn.is_balanced()


################################################################################
# Negative controls -- what must STILL be refused
################################################################################

def test_electron_gain_without_charge_change_is_still_refused():
    """THE control.  ``Ar + e- => Ar + e- + e-`` invents an electron out of
    nothing: the heavy census matches, but charge goes -1 -> -2.

    Removing the electron from the element census must not remove it from the
    check.  This is refused before the change (on the census) and after it (on
    charge), and it is what a reviewer should try first.
    """
    rxn = _rxn(['Ar', 'e'], ['Ar', 'e', 'e'])
    assert _charges(rxn) == (-1, -2)
    assert not rxn.is_balanced()


def test_electron_loss_without_charge_change_is_still_refused():
    """The same control in the other direction: ``Ar + e- => Ar`` destroys an
    electron, charge -1 -> 0."""
    rxn = _rxn(['Ar', 'e'], ['Ar'])
    assert _charges(rxn) == (-1, 0)
    assert not rxn.is_balanced()


def test_ionisation_with_the_wrong_number_of_electrons_is_still_refused():
    """``Ar + e- => Ar+ + e- + e- + e-`` ionises once but emits two electrons:
    charge -1 -> -2.  A single-electron-off error in an otherwise plausible
    ionisation is still caught."""
    rxn = _rxn(['Ar', 'e'], ['Arplus', 'e', 'e', 'e'])
    assert _charges(rxn) == (-1, -2)
    assert not rxn.is_balanced()


def test_charge_non_conserving_reaction_without_electrons_is_still_refused():
    """``N2+ + N2 => N2 + N + N``, carried verbatim in ``PlasmaAir`` and
    charge-broken at source: +1 -> 0, with no free electron anywhere in it.

    It is the sole ``PlasmaAir`` entry this change does not repair, and it is
    green both before and after.  A change that turns it green has switched the
    check off rather than fixed it.
    """
    rxn = _rxn(['N2plus', 'N2'], ['N2', 'N', 'N'])
    assert _charges(rxn) == (1, 0)
    assert not rxn.is_balanced()


def test_heavy_element_mismatch_is_still_refused():
    """The census still runs for every element that is not ``e``: ``Ar + e- =>
    Ar + Ar+ + e- + e-`` conserves charge (-1 = -1) but invents an argon atom."""
    rxn = _rxn(['Ar', 'e'], ['Ar', 'Arplus', 'e', 'e'])
    assert _charges(rxn) == (-1, -1)
    assert not rxn.is_balanced()


def test_scalar_electrons_are_still_folded_into_the_charge_comparison():
    """``Ar => Ar+`` with ``electrons=+1`` balances through the scalar path, and
    the same equation with ``electrons=0`` does not.  The census skip does not
    bypass ``Reaction.electrons``."""
    assert _rxn(['Ar'], ['Arplus'], electrons=1).is_balanced()
    assert not _rxn(['Ar'], ['Arplus'], electrons=0).is_balanced()


################################################################################
# Why the widening is narrow
################################################################################

@pytest.mark.parametrize('reactants,products,delta_e,delta_heavy_charge', [
    (['Ar', 'e'], ['Arplus', 'e', 'e'], +1, +1),
    (['Ar'], ['Arplus', 'e'], +1, +1),
    (['Arplus', 'e'], ['Ar'], -1, -1),
    (['Arplus', 'e', 'e'], ['Ar', 'e'], -1, -1),
])
def test_charge_conserving_electron_change_implies_matching_ionisation(
        reactants, products, delta_e, delta_heavy_charge):
    """Every reaction this change newly admits is an ionisation or attachment by
    exactly the number of electrons it gains or loses.

    The free electron is the only participant whose element is ``e`` and it
    carries charge -1 exactly, so conserving total charge while the electron
    count moves by ``delta_e`` forces the heavy species' charge to move by
    ``+delta_e``.  There is no room for a reaction with a genuine
    electron-accounting error that nonetheless conserves charge -- which is the
    question ``is_balanced``'s reviewers should be asking, answered by
    construction rather than by enumeration.
    """
    rxn = _rxn(reactants, products)

    def counts(species_list):
        n_e = sum(1 for spc in species_list if spc.molecule[0].atoms[0].element.symbol == 'e')
        heavy_charge = sum(atom.charge for spc in species_list
                           for atom in spc.molecule[0].atoms
                           if atom.element.symbol != 'e')
        return n_e, heavy_charge

    r_e, r_q = counts(rxn.reactants)
    p_e, p_q = counts(rxn.products)
    assert p_e - r_e == delta_e
    assert p_q - r_q == delta_heavy_charge
    assert delta_heavy_charge == delta_e
    assert rxn.is_balanced()
