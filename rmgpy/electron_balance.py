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
    'get_plasma_rate_order',
    'potential_dependence_is_inert',
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
    positive means they are produced (product side). ``Reaction.electrons`` is
    signed relative to the reaction object's current reactant/product
    orientation, so reversing the object negates it; this helper reads that
    current-orientation sign. (``KineticsFamily.electrons`` is a different thing --
    the family-forward declaration.)

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


def potential_dependence_is_inert(kinetics):
    """
    Return ``True`` when a charge-transfer rate law's potential dependence is
    identically absent, so that writing its ``V = V0`` rate ``(A, n, Ea)`` is
    exact at every potential rather than only at the reference one.

    Both :class:`SurfaceChargeTransfer` and :class:`ArrheniusChargeTransfer`
    evaluate::

        k(T, V) = A * (T/T0)^n * exp(-Ea_eff / (R*T))
        Ea_eff  = Ea - alpha * electrons * F * (V - V0)

    so ``Ea_eff`` is free of ``V`` exactly when ``alpha * electrons == 0``.

    The second condition is less obvious and just as load-bearing.
    ``get_activation_energy_from_potential`` clamps a negative ``Ea_eff`` up to
    zero, and ``get_rate_coefficient`` only calls it on the ``V != V0`` branch --
    so a rate with ``Ea < 0`` still jumps as soon as ``V`` leaves ``V0``, even
    when ``alpha * electrons == 0``. Measured at ``alpha=0``, ``electrons=-1``,
    ``Ea=-5 kJ/mol``: ``k(V0)=3.33e+06`` against ``k(V!=V0)=1.00e+06``.

    Anything else is a live potential dependence with no exact reduction, and the
    writers refuse it rather than emit the reference-potential number.
    """
    alpha = getattr(kinetics, 'alpha', None)
    electrons = getattr(kinetics, 'electrons', None)
    Ea = getattr(kinetics, 'Ea', None)
    if alpha is None or electrons is None or Ea is None:
        return False
    return (alpha.value_si * electrons.value_si == 0.0) and Ea.value_si >= 0.0


#: Rate-coefficient units to the reaction order they imply. Same mapping the
#: kinetics classes use in their own ``to_cantera_kinetics``.
_ORDER_BY_RATE_UNITS = {
    '1/s': 1,
    's^-1': 1,
    'm^3/(mol*s)': 2,
    'cm^3/(mol*s)': 2,
    'm^3/(molecule*s)': 2,
    'cm^3/(molecule*s)': 2,
    'm^6/(mol^2*s)': 3,
    'cm^6/(mol^2*s)': 3,
    'm^6/(molecule^2*s)': 3,
    'cm^6/(molecule^2*s)': 3,
}


def get_plasma_rate_order(kinetics):
    """
    Return the reaction order implied by a plasma rate coefficient, or ``None``
    when the kinetics is not a plasma type or its units are not recognised.

    :class:`ElectronCollisionPlasma` stores a cross-section rather than an
    A-factor; its rate coefficient is the Maxwellian average ``<sigma*v>``, which
    is bimolecular by construction, so it is always order 2.
    """
    from rmgpy.kinetics import (
        TwoTemperaturePlasma, ElectronCollisionPlasma,
        BadnellRRArrhenius, VoronovEIArrhenius,
    )

    if isinstance(kinetics, ElectronCollisionPlasma):
        return 2

    if isinstance(kinetics, (TwoTemperaturePlasma, BadnellRRArrhenius, VoronovEIArrhenius)):
        A = getattr(kinetics, 'A', None)
        return _ORDER_BY_RATE_UNITS.get(getattr(A, 'units', None))

    return None


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
    if kinetics is None:
        return

    if any(spc.is_electron() for spc in reactants):
        return

    if getattr(kinetics, 'uses_electron_density', False):
        # BadnellRRArrhenius and VoronovEIArrhenius say so themselves.
        raise MechanismWriterError(
            'Reaction {0!s} has kinetics {1} whose rate is proportional to the electron '
            'density, but the exported equation "{2}" has no electron among its reactants, '
            'so a solver would evaluate it at the wrong reaction order. Put the consumed '
            'electron in reaction.reactants and count only surplus produced electrons in '
            'reaction.electrons.'.format(reaction, type(kinetics).__name__, equation)
        )

    # TwoTemperaturePlasma and ElectronCollisionPlasma do not carry
    # uses_electron_density at all, so asking the flag is useless for exactly the
    # two classes most likely to trip this. Ask the rate coefficient instead: its
    # dimensionality fixes the reaction order, and a plasma rate coefficient that
    # is one order higher than the exported reactant side is missing its electron.
    required_order = get_plasma_rate_order(kinetics)
    if required_order is None or required_order == len(reactants):
        return

    raise MechanismWriterError(
        'Reaction {0!s} has kinetics {1} whose rate coefficient is of order {2:d}, but the '
        'exported equation "{3}" has {4:d} reactant(s) and no electron among them, so a solver '
        'would evaluate it at the wrong reaction order. An electron-impact reaction must carry '
        'its electron in reaction.reactants (or in reaction.electrons as a negative count).'
        .format(reaction, type(kinetics).__name__, required_order, equation, len(reactants))
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
