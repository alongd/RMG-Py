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
Shared helpers for getting the electron stoichiometry of a charged reaction into
an exported mechanism, used by both the Chemkin and the Cantera writers.

RMG carries the electron stoichiometry of a charged reaction in the scalar
``Reaction.electrons``, not in ``Reaction.reactants``/``Reaction.products``.
Both writers build their equation strings from the reactant and product lists
alone, so without these helpers ``reaction.electrons`` is never serialized and
every electron-consuming or electron-producing reaction is exported
element-unbalanced in the ``E`` pseudo-element.

Both writers must behave identically here, so the logic lives in one place
rather than being written twice -- once in Python and once in Cython.
"""

from rmgpy.exceptions import MechanismWriterError

__all__ = [
    'get_electron_species',
    'expand_electrons',
    'get_species_electron_count',
    'check_electron_balance',
    'check_electron_reactant_order',
]


def get_electron_species(species_list):
    """
    Return the electron :class:`Species` in ``species_list``, or ``None`` if the
    mechanism has no electron species.
    """
    for spc in species_list or []:
        try:
            if spc.is_electron():
                return spc
        except (AttributeError, IndexError):
            continue
    return None


def expand_electrons(reaction, species_list):
    """
    Return ``(reactants, products)``: copies of ``reaction.reactants`` and
    ``reaction.products`` with ``reaction.electrons`` folded in as an explicit
    electron species, so that the exported equation carries the electron the way
    a solver needs to see it.

    Sign convention matches :meth:`rmgpy.reaction.Reaction.is_balanced`: negative
    ``electrons`` means electrons are consumed (they belong on the reactant side),
    positive means they are produced (product side).

    Raises :class:`MechanismWriterError` if the reaction needs an electron but the
    mechanism does not define an electron species -- exporting the equation without
    it would be exactly the silent corruption these helpers exist to prevent.
    """
    reactants = list(reaction.reactants)
    products = list(reaction.products)

    electrons = getattr(reaction, 'electrons', 0) or 0
    if electrons == 0:
        return reactants, products

    electron = get_electron_species(species_list)
    if electron is None:
        raise MechanismWriterError(
            'Reaction {0!s} has electrons={1:d} but the mechanism defines no electron '
            'species, so the electron cannot be written into the exported equation. '
            'Add the electron to the species list before exporting.'.format(reaction, electrons)
        )

    if electrons < 0:
        reactants.extend([electron] * (-electrons))
    else:
        products.extend([electron] * electrons)

    return reactants, products


def get_species_electron_count(species):
    """
    Return the number of ``E`` pseudo-element units a species contributes, using
    the same rule both writers use to build the species composition block: ``E``
    is minus the net charge for a charged species, otherwise the count of
    electron atoms in the structure.

    Cantera's convention: negatively charged ions have ``E > 0``, positively
    charged ions have ``E < 0``.
    """
    molecule = species.molecule[0]
    charge = molecule.get_net_charge()
    if charge != 0:
        return -charge
    return sum(1 for atom in molecule.atoms if atom.element.chemkin_name == 'E')


def check_electron_reactant_order(reaction, reactants, equation):
    """
    Raise :class:`MechanismWriterError` if the reaction's rate law is proportional
    to the electron density but the exported reactant side has no electron on it.

    An electron-impact ionization written as ``A => A+ + e`` (net electrons = +1)
    balances in ``E`` and looks fine, but a solver reads it as first order in A
    and never multiplies by the electron concentration -- so the rate is wrong by
    a factor of the electron density while the file looks well formed. The
    correct RMG representation puts the consumed electron in
    ``reaction.reactants`` and counts only the surplus produced electrons in
    ``reaction.electrons`` (``A + e => A+ + 2 e`` is ``reactants=[A, e]`` with
    ``electrons=2``).
    """
    kinetics = getattr(reaction, 'kinetics', None)
    if not getattr(kinetics, 'uses_electron_density', False):
        return
    if any(spc.is_electron() for spc in reactants):
        return
    raise MechanismWriterError(
        'Reaction {0!s} has kinetics {1} whose rate is proportional to the electron '
        'density, but the exported equation "{2}" has no electron among its reactants, '
        'so a solver would evaluate it at the wrong reaction order. Put the consumed '
        'electron in reaction.reactants and count only surplus produced electrons in '
        'reaction.electrons.'.format(reaction, type(kinetics).__name__, equation)
    )


def check_electron_balance(reaction, reactants, products, equation):
    """
    Raise :class:`MechanismWriterError` unless the ``E`` pseudo-element balances
    across ``equation``, whose sides are the already-electron-expanded
    ``reactants`` and ``products`` lists.

    This is checked on what is about to be written, not on the RMG objects, so it
    catches the case where ``Reaction.is_balanced`` and the writer disagree -- and
    only the writer reaches the solver.
    """
    reactant_e = sum(get_species_electron_count(spc) for spc in reactants)
    product_e = sum(get_species_electron_count(spc) for spc in products)
    if reactant_e != product_e:
        raise MechanismWriterError(
            'Reaction {0!s} does not balance in the E pseudo-element and would be '
            'exported wrong: the equation "{1}" has E={2:d} on the left and E={3:d} on '
            'the right (reaction.electrons={4:d}).'.format(
                reaction, equation, reactant_e, product_e, getattr(reaction, 'electrons', 0) or 0)
        )
