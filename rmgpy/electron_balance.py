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

That "one place" is also where the export path meets the *other* consumer of the
same conversion, :class:`rmgpy.solver.plasma.PlasmaReactor`, and I-126 is the
ticket that joined them. The reactor converts through
:func:`rmgpy.electron_placement.resolve_electron_placement`, driven by the
two-sided ``FAMILY_ELECTRON_PLACEMENT`` declaration; the export path converted
through :func:`expand_electrons`, driven by the net scalar alone. Two converters,
one representation, and they disagreed on any shape the net scalar cannot
express -- so a mechanism the reactor evaluated correctly could not be written
out at all. :func:`expand_electrons` now reads that same declaration wherever the
reaction's owner carries one, and falls back to the net rule only where no owner
has spoken. The declaration is stated once, in ``electron_placement``, and both
halves of the codebase read it.
"""

from rmgpy.exceptions import MechanismWriterError

__all__ = [
    'get_electron_species',
    'get_placement_owner',
    'get_placement_declaration',
    'get_electron_placement_counts',
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


def get_placement_owner(reaction):
    """
    Return the label whose :data:`FAMILY_ELECTRON_PLACEMENT` declaration governs
    ``reaction``, or ``None`` when no attribution the reaction carries names a
    declared owner.

    Two attributions are consulted, in a fixed precedence, and both travel WITH
    the reaction object rather than being read off whatever container currently
    holds it:

    * ``Reaction.family`` -- the current attribution: the family label for a
      :class:`TemplateReaction`, the library label for a
      :class:`LibraryReaction` (``LibraryReaction.__init__`` sets
      ``family = library``), or the fixed container label (``seed`` /
      ``restart``) for a reaction reloaded from an auto-generated seed
      mechanism, because ``KineticsLibrary.get_library_reactions`` overwrites
      ``family`` with the label of the library it loaded from.
    * ``Reaction.library`` -- the reaction's ORIGINAL owner, for exactly the
      reloaded-seed case above: the seed writer emits ``Originally from
      reaction library: <label>`` into each entry's ``longDesc``
      (``rmgpy/rmg/main.py``, ``make_seed_mech``), and
      ``get_library_reactions`` parses that line back into
      ``LibraryReaction.library`` for auto-generated libraries. This is the
      I-148 repair: before it, only ``family`` was consulted, so every plasma
      reaction that round-tripped through a seed lost its declaration to the
      renamed container and silently fell back to the net-derived rule -- the
      ionisation channel's (1, 2) collapsed to (0, 1), the restarted core
      carried the channel twice, and the first mechanism save was refused.

    ``family`` wins when both name a declared owner: the current attribution
    governs, and the preserved provenance is consulted only where the current
    attribution is silent. This rule is monotonic over the pre-I-148 behaviour
    -- every reaction that resolved before resolves to the same owner, and the
    only reactions that newly resolve are those whose ``family`` names no
    declaration at all (the seed's fixed labels are never declared; declaring
    them is the one-line non-fix this module refuses, since it would hand one
    placement to every reaction that ever passed through a seed).

    A plain :class:`Reaction` -- which is what the reactor's placement view is
    -- has neither attribute, so a view can never be re-expanded through a
    declaration.

    Imported lazily, like :func:`get_placement_declaration` and for the same
    reason.
    """
    from rmgpy.electron_placement import FAMILY_ELECTRON_PLACEMENT

    family = getattr(reaction, 'family', None)
    if family and family in FAMILY_ELECTRON_PLACEMENT:
        return family
    library = getattr(reaction, 'library', None)
    if library and library in FAMILY_ELECTRON_PLACEMENT:
        return library
    return None


def get_placement_declaration(reaction):
    """
    Return the ``(reactant_count, product_count)`` electron-placement declaration
    the reaction's owner carries, or ``None`` if no attribution the reaction
    carries names an owner declared in
    :data:`rmgpy.electron_placement.FAMILY_ELECTRON_PLACEMENT`.

    The owner is resolved by :func:`get_placement_owner`, which reads the
    reaction's own attributions -- ``family`` first, then the ``library``
    provenance a seed round trip preserves -- never the label of the container
    the reaction currently sits in.

    ``None`` means "this owner has made no statement about placement", which is
    the case for every reaction in RMG except the handful of declared plasma
    owners, and which leaves :func:`expand_electrons` on its net-derived rule.
    It never means "the declaration was unusable": a malformed declaration
    raises, because silently falling back to the very rule the declaration
    exists to override is how a wrong equation gets written while the export
    reports success.

    Imported lazily: :mod:`rmgpy.electron_placement` imports this module at
    module scope, and a name lookup inside the function keeps that a one-way
    dependency.
    """
    from rmgpy.electron_placement import FAMILY_ELECTRON_PLACEMENT

    owner = get_placement_owner(reaction)
    if owner is None:
        return None
    declaration = FAMILY_ELECTRON_PLACEMENT[owner]
    if (not isinstance(declaration, tuple) or len(declaration) != 2
            or not all(isinstance(n, int) and not isinstance(n, bool) and n >= 0
                       for n in declaration)
            or declaration == (0, 0)):
        raise MechanismWriterError(
            'Owner {0!r} carries the malformed electron-placement declaration {1!r} '
            '(reaction {2!s}). A declaration must be a (reactant_count, product_count) '
            'pair of non-negative integers placing at least one electron; refusing to '
            'export the reaction rather than fall back to the net-derived placement the '
            'declaration exists to override.'.format(owner, declaration, reaction)
        )
    return declaration


def get_electron_placement_counts(reaction):
    """
    Return ``(reactant_count, product_count)``: how many free electrons stand on
    each side of ``reaction`` **as it is currently oriented**.

    This is the electron half of reaction IDENTITY. Two reactions are the same
    reaction when they have the same participants on each side, and the electron
    is a participant like any other -- it is merely stored as metadata rather
    than in the participant lists. This function restores it to a comparable
    form so that identity predicates can include it.

    **Why this is not ``Reaction.electrons``, which is the trap.** The net scalar
    is ONE number and identity needs TWO, for exactly the reason
    :func:`expand_electrons` needs two: a net count cannot say how many electrons
    were incident, only how many appeared or vanished. The consequence for
    identity is sharper than the consequence for export, and it is what
    ``docs/i134-duplicate-electrons/`` exists to record. Take the two shipped
    lithium plasma channels:

        electron-impact ionisation   Li  + e-  =>  Li+ + 2 e-    placement (1, 2)
        radiative recombination      Li+ + e-  =>  Li  + hv      placement (1, 0)

    Their heavy species are exact mirrors, so every reference comparison matches
    them in the reverse direction. Their NET electron counts are ``+1`` and
    ``-1`` -- exactly equal and opposite, which is precisely the relation a
    reverse pair has. **So a net-based test declares them the same reaction, and
    they are not.** The reverse of the ionisation is three-body recombination,
    ``Li+ + 2 e- => Li + e-``, placement ``(2, 1)``, third order; the radiative
    channel is ``(1, 0)``, second order, with a photon carrying off the energy.
    Different molecularity, different rate law, different temperature dependence,
    opposite roles in the charge balance. Only the placement pair separates them:
    ``(1, 2)`` reversed is ``(2, 1)``, which is not ``(1, 0)``.

    This is the same distinction ``FAMILY_ELECTRON_PLACEMENT`` was widened for in
    I-113 and I-126 -- incident order is DECLARED, never derived -- turning out to
    be what identity needs as well. Anyone reaching for ``reaction.electrons``
    here should read the paragraph above first: it is not equivalent, and the
    failure it produces is silent.

    **The contract for reactions with no declaration**, which is almost every
    reaction in RMG. When the owner has made no statement, the counts come from
    the same net-derived rule :func:`expand_electrons` falls back to: negative
    ``electrons`` means they are consumed, so ``(-electrons, 0)``; positive means
    produced, so ``(0, electrons)``; zero means ``(0, 0)``. Two consequences of
    that rule are worth stating precisely, because they bound what consulting
    this function can change:

    * it reduces to a comparison of the net counts. ``(r1, p1) == (r2, p2)``
      holds exactly when ``electrons1 == electrons2``, since one count is always
      zero and the other is ``|electrons|`` on the side the sign names; and
      ``(r1, p1) == (p2, r2)`` holds exactly when ``electrons1 == -electrons2``,
      by the same argument with the sides swapped. Those are precisely the two
      comparisons :meth:`rmgpy.reaction.Reaction.is_isomorphic` already made on
      the net scalar, so **that** predicate's verdict for an undeclared owner is
      unchanged bit for bit.
    * ``electrons = 0`` gives ``(0, 0)``, so any two reactions that are both
      electrically neutral compare equal on electrons whatever else is true of
      them. That is what bounds the change for
      :func:`rmgpy.rmg.model.are_identical_species_references`, which compared no
      electron information at all before: for neutral chemistry -- all of RMG
      outside the charged families and the plasma libraries -- its verdict is
      also unchanged. For a CHARGED reaction it is a genuine tightening, which is
      the point of the repair.

    A declared owner whose declaration is ONE-SIDED -- ``(n, 0)`` or ``(0, n)``
    -- places exactly where the net rule would, so it moves nothing either. Only
    a two-sided declaration can change an answer that the net rule got, and
    ``PlasmaElectronImpactIonization``'s ``(1, 2)`` is the only one shipped.

    **Orientation.** A declaration is stated in its owner's forward orientation.
    A reaction object may be stored either way round, and ``Reaction.electrons``
    is signed to that stored orientation, so the two together determine it: the
    stored orientation is forward when ``product_count - reactant_count`` equals
    ``electrons``, and reversed when ``reactant_count - product_count`` does. The
    reversed case returns the declaration swapped.

    **When the reaction and its owner disagree**, matching neither orientation,
    this falls through to the net rule rather than raising.
    :func:`expand_electrons` raises there, and should: it is an export boundary,
    and writing a wrong equation is unrecoverable. This function is consulted on
    every duplicate comparison of every run, where raising would convert a data
    inconsistency into a traceback from an unrelated-looking place, and where the
    export boundary is still downstream and still refuses. Falling back is
    reading the reaction's own account of itself when its owner's account cannot
    be applied. No shipped path reaches this branch: ``_create_reaction`` sets
    ``electrons`` from the family declaration itself, and library entries are
    balance-checked at load.

    ``reaction`` is never mutated, and explicit electrons already sitting in the
    participant lists are not counted -- a reaction in that form (the reactor's
    placement view, or a mechanism read back from a written file) carries
    ``electrons = 0`` and has its electrons compared as ordinary participants by
    the list comparison itself.
    """
    electrons = getattr(reaction, 'electrons', 0) or 0

    declaration = get_placement_declaration(reaction)
    if declaration is not None:
        reactant_count, product_count = declaration
        if product_count - reactant_count == electrons:
            return reactant_count, product_count
        if reactant_count - product_count == electrons:
            return product_count, reactant_count

    if electrons < 0:
        return -electrons, 0
    return 0, electrons


def expand_electrons(reaction, species_list):
    """
    Return ``(reactants, products)``: copies of ``reaction.reactants`` and
    ``reaction.products`` with ``reaction.electrons`` folded in as an explicit
    electron species, so that the exported equation carries the electron the way
    a solver needs to see it.

    There are two placement rules, and which one applies is decided by the
    reaction's OWNER, never by the writer.

    **The declaration, when the owner has one.** If the owner (family or kinetics
    library) appears in ``FAMILY_ELECTRON_PLACEMENT``, its
    ``(reactant_count, product_count)`` pair is authoritative and the electron is
    placed on both sides in the numbers it names. This is the same declaration,
    read from the same table, that
    :func:`rmgpy.electron_placement.resolve_electron_placement` gives the plasma
    reactor -- which is the point: the reactor and the writers then cannot
    disagree about which representation a reaction is in, because there is only
    one statement of it. ``reaction.electrons`` is validated against the
    declaration's net change (``product_count - reactant_count``) and read for
    nothing else.

    **The net-derived rule, when it has not.** Negative ``electrons`` means
    electrons are consumed (they belong on the reactant side), positive means
    they are produced (product side), matching
    :meth:`rmgpy.reaction.Reaction.is_balanced`. ``Reaction.electrons`` is signed
    relative to the reaction object's current reactant/product orientation, so
    reversing the object negates it; this helper reads that current-orientation
    sign. (``KineticsFamily.electrons`` is a different thing -- the
    family-forward declaration.)

    Why the net rule is not enough on its own, which is what I-126 fixed: a net
    count is ONE number and placement needs TWO. ``electrons = +1`` is equally
    consistent with ``Li + e- => Li+ + 2 e-`` (order 2, the
    ``cm^3/(molecule*s)`` Voronov coefficient's own dimensionality) and with
    ``Li => Li+ + e-`` (order 1), and the net rule always produces the second.
    Both balance in ``E``, so the resulting file looks well formed while the rate
    is wrong by a factor of the electron density -- which is exactly what
    :func:`check_electron_reactant_order` refused to write, correctly. The net
    rule stays for every owner that has made no declaration, where it is not
    merely adequate but provably identical to the declared one: an owner whose
    declaration is one-sided (``(n, 0)`` or ``(0, n)``) places exactly where the
    net count would.

    Raises :class:`MechanismWriterError` if the reaction needs an electron but the
    mechanism does not define an electron species -- exporting the equation without
    it would be exactly the silent corruption these helpers exist to prevent --
    and if a declared owner's reaction carries a net count its declaration does
    not predict.
    """
    reactants = list(reaction.reactants)
    products = list(reaction.products)

    electrons = getattr(reaction, 'electrons', 0) or 0
    if electrons == 0:
        # Nothing to fold in. A reaction already in explicit-electron form (the
        # reactor's placement view, or a mechanism read back from a written
        # file) arrives here and is returned untouched, which is what makes this
        # helper idempotent over the write/read round trip.
        return reactants, products

    electron = get_electron_species(species_list)
    if electron is None:
        raise MechanismWriterError(
            'Reaction {0!s} has electrons={1:d} but the mechanism defines no electron '
            'species, so the electron cannot be written into the exported equation. '
            'Add the electron to the species list before exporting.'.format(reaction, electrons)
        )

    declaration = get_placement_declaration(reaction)
    if declaration is not None:
        reactant_count, product_count = declaration
        if product_count - reactant_count != electrons:
            raise MechanismWriterError(
                'Reaction {0!s} carries electrons={1:+d}, but its owner {2!r} declares the '
                'electron placement {3!r} -- net {4:+d} ({5:d} electron(s) on the reactant '
                'side, {6:d} on the product side). The reaction and its owner describe '
                'different chemistry; refusing to export it rather than prefer either '
                'source.'.format(reaction, electrons, getattr(reaction, 'family', None),
                                 declaration, product_count - reactant_count,
                                 reactant_count, product_count)
            )
        reactants.extend([electron] * reactant_count)
        products.extend([electron] * product_count)
        return reactants, products

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


def is_isomorphic_same_charge(species, other, strict=False):
    """
    Return ``True`` iff ``species`` and ``other`` are isomorphic under
    ``is_isomorphic(..., strict=strict)`` AND carry the same net charge.

    ``strict=False`` isomorphism is charge-blind (``Atom.equivalent`` with
    ``strict=False`` compares only the element), so a single O atom of charge -1
    is isomorphic to one of charge +1. Several lookups match with ``strict=False``
    on purpose, so that resonance structures of a species match each other -- but
    on their own they then treat an anion and a cation with the same heavy-atom
    skeleton as the same species. Net charge is conserved across resonance
    structures, so requiring equal net charge never rejects a genuine resonance
    match; it only stops opposite-charge species from being confused.

    ``species`` must be an :class:`rmgpy.species.Species` (it owns the
    ``is_isomorphic`` call); ``other`` may be a :class:`~rmgpy.species.Species` or
    a :class:`~rmgpy.molecule.molecule.Molecule`. Both expose ``get_net_charge``.
    """
    if species.get_net_charge() != other.get_net_charge():
        return False
    return species.is_isomorphic(other, strict=strict)


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
