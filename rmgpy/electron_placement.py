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
The electron-representation boundary: family-declared electron placement.

RMG's database layer carries the electron of an electron-consuming reaction as
the scalar ``Reaction.electrons`` (net stoichiometry metadata), not as an
explicit participant. The plasma reactor deliberately refuses that
representation: a net count cannot express incident-electron order, so
:class:`rmgpy.solver.plasma.PlasmaReactor` rejects any reaction with nonzero
``electrons`` and instead counts incident electrons from explicit reactant
placement (``Species.is_electron()``).

This module is the authorised upstream converter between the two
representations. :func:`resolve_electron_placement` produces a non-mutating,
reactor-facing VIEW of a reaction: copied reactant/product lists with the
canonical electron species inserted on each side in the numbers the reaction's
*family* declares, and ``electrons = 0`` in the view. The canonical database reaction is
never modified — it keeps its metadata count and its original participant
lists, so the database representation (and the training-set label convention
built on it) is untouched.

The order source is the family-level placement declaration in
:data:`FAMILY_ELECTRON_PLACEMENT`, established fresh each time the resolver
runs. It is NEVER the net scalar: ``electrons`` is only *validated against* the
declaration, and the placement itself is performed by the local list-append
primitive :func:`_place_declared_electrons`, which appends exactly the declared
counts and knows nothing about the net scalar, with the result checked against
the declaration afterwards.

The declaration is two-sided — ``(reactant_count, product_count)`` — because the
electron's contribution to reaction ORDER and its NET stoichiometric change are
two different numbers, and only attachment makes them coincide. ``A + e- -> A-``
is one electron in, none out: that single number is both the whole net change
and the whole of the electron's contribution to order, which is why a one-sided
``(side, count)`` declaration served for as long as attachment was the only
shape. Ionisation separates them: ``Li + e- -> Li+ + 2e-`` has incident order 1
— which fixes the rate as second order overall and is what makes the
``cm^3/(molecule*s)`` Voronov coefficient dimensionally meaningful — and net
production +1. A bare ``electrons = +1`` cannot tell that apart from
``Li -> Li+ + e-``, which is first order and would need ``s^-1``; get it wrong
and the rate is off by a factor of the electron density while the file still
looks well formed. So incident order is ``reactant_count``, *declared*, and net
is ``product_count - reactant_count``, *validated*.

Scope: the owners named in :data:`FAMILY_ELECTRON_PLACEMENT` — families and
kinetics libraries alike, since ``LibraryReaction`` sets ``family = library`` and
the resolver reads ``Reaction.family`` — in the shape each declares. Every other
owner, reaction shape, or kinetics form fails by
name with :class:`rmgpy.exceptions.ElectronPlacementError` — there is no silent
fallback to net-derived inference, and no general mechanism to inherit
accidentally.

Direction and reversibility are deliberately NOT part of that scope, since
I-088. The declaration is about the family-forward molecular orientation, and the
engine stores every generated reaction in that orientation whichever direction it
was generated in, so ``is_forward`` says nothing about which side the electron
belongs on; the resolver verifies the side structurally instead (step 10, the
``E`` pseudo-element balance of the finished view). Reversibility is a
reverse-RATE question, owned and enforced by the consumer —
:class:`rmgpy.solver.plasma.PlasmaReactor` refuses a reversible
electron-containing reaction by name — and the view carries ``reversible``
through unchanged so that guard still fires.
"""

import logging

from rmgpy.electron_balance import get_plasma_rate_order, get_species_electron_count
from rmgpy.exceptions import ElectronPlacementError
from rmgpy.reaction import Reaction

__all__ = [
    'FAMILY_ELECTRON_PLACEMENT',
    'RATE_ORDER_AGREES',
    'RATE_ORDER_UNRESOLVABLE',
    'resolve_electron_placement',
]

#: Family-level electron placement declarations, keyed by family label.
#: Each value is ``(reactant_count, product_count)``: how many electrons the
#: family places on each side of the reaction, in its FAMILY-FORWARD
#: orientation. This mapping — not ``Reaction.electrons`` — is the order source
#: for :func:`resolve_electron_placement`. A family absent from this mapping has
#: NO placement declaration and resolves to a named failure, never to a
#: net-derived guess.
#:
#: The two numbers are independent and neither is derivable from the other.
#: ``reactant_count`` is the INCIDENT ORDER — the electron's contribution to the
#: reaction order, which is what the rate coefficient's dimensionality has to
#: match. ``product_count - reactant_count`` is the NET change, which is the only
#: thing ``Reaction.electrons`` is ever compared against. Attachment is the
#: degenerate case where one number does both jobs; ionisation is not (see the
#: module docstring).
#:
#: Three of the four declared owners place a single electron on the reactant side
#: and none on the product side, for three unrelated chemistries — which is
#: exactly why ``(1, 0)`` has to be declared three times rather than inferred once:
#: ``Plasma_Electron_Attachment`` is non-dissociative attachment
#: (``A + e- -> A-``); ``Cation_R_Recombination``
#: is cation-radical recombination (``Li+ + R. + e- -> R-Li``), which the plasma
#: decks reach in the family's REVERSE generation direction — RMG has only the
#: neutral ``R-Li``, so it matches the product template and reconstructs
#: ``Li+ + R.``. The declaration is about the family's forward molecular
#: orientation, which is the orientation the engine stores in BOTH generation
#: directions (see step 6 below), so one entry covers both.
#:
#: ``PlasmaRadiativeRecombination`` is the third, added by I-119: ``A+ + e- ->
#: A + hv``, the Badnell (2006) fits, the electron SINK that pairs with the
#: ionisation source below. It shares attachment's numbers because it shares
#: attachment's SHAPE — the electron is captured and does not come out again, so
#: the one electron is both the whole incident order and the whole net change.
#: Sharing them buys it nothing: an owner absent from this mapping resolves to a
#: named failure even when another owner has already declared its exact shape,
#: which ``test_plasma_radiative_recombination.py`` in RMG-database pins. The
#: photon is not represented — it carries neither charge nor mass, so no balance
#: check misses it.
#:
#: One thing that entry is deliberately NOT: the reverse of the ionisation entry
#: below. Radiative recombination is a second-order forward channel with a photon;
#: the true inverse of electron-impact ionisation is three-body recombination
#: ``A+ + 2 e- -> A + e-``, which would declare ``(2, 1)``. That declaration is
#: absent because the channel cannot be stored at all — ``TwoTemperaturePlasma``,
#: the only Te-aware third-order rate law, carries no ``electrons`` field, so such
#: an entry fails ``KineticsLibrary.load``'s balance check before any question of
#: data arises. Measured, with what would be needed to lift it, in
#: ``docs/i119-recombination-loss.md`` in RMG-database — which also carries the
#: coverage arithmetic that narrows 318 Badnell stages to the single usable Li+.
#:
#: The fourth, ``PlasmaElectronImpactIonization``, declares ``(1, 2)`` — one
#: electron incident, two liberated — and is the only declaration whose two
#: numbers differ, which is the case I-113 widened the schema for. It is a
#: kinetics LIBRARY label, not a family label, and that is deliberate rather than
#: accidental: what this resolver reads is ``Reaction.family``, and
#: ``LibraryReaction.__init__`` sets ``self.family = library``
#: (``rmgpy/data/kinetics/library.py``), so a library owns its placement on
#: exactly the same terms a family does. The mapping's name is simply older than
#: the question. A family was rejected on measurement rather than taste — RMG's
#: rate-rule averaging cannot combine two ``VoronovEIArrhenius`` fits at all, so
#: a Voronov-trained family is capped at one training entry per node and buys no
#: generativity; the argument is ``docs/i114-ionisation-owner.md`` in
#: RMG-database.
#:
#: The registry stays a closed, hand-maintained list: an entry appearing here
#: means someone decided it should.
FAMILY_ELECTRON_PLACEMENT = {
    'Plasma_Electron_Attachment': (1, 0),
    'Cation_R_Recombination': (1, 0),
    'PlasmaRadiativeRecombination': (1, 0),
    'PlasmaElectronImpactIonization': (1, 2),
}

#: Outcomes of the rate-order cross-check (step 11). Disagreement is not one of
#: them — it raises. These two are recorded on the view rather than returned,
#: and they are deliberately DISTINCT: ``get_plasma_rate_order`` answers ``None``
#: for any rate law whose units are not in its table, and a silent ``None`` must
#: never read as agreement.
RATE_ORDER_AGREES = 'rate-order cross-check: agrees'
RATE_ORDER_UNRESOLVABLE = 'rate-order cross-check: order not resolvable'

#: Kinetics classes whose own ``electrons`` field is documented to be the SAME
#: signed net quantity as ``Reaction.electrons`` — "the net change in free
#: electron number, signed to the reaction as written" — and which can therefore
#: be validated against the declaration exactly as the reaction scalar is,
#: rather than refused as a second placement source.
#:
#: Every OTHER kinetics-level ``electrons`` stays fatal. The charge-transfer rate
#: laws carry one meaning "the number of electrons transferred", which is a
#: different quantity that no declaration in this table speaks to; reconciling it
#: would be guessing. Kept as a closed list for the same reason
#: :data:`FAMILY_ELECTRON_PLACEMENT` is.
_NET_ELECTRON_KINETICS_CLASSES = ('BadnellRRArrhenius', 'VoronovEIArrhenius')


def _is_electron(spc):
    """
    ``True`` if ``spc`` (a :class:`Species` or :class:`Molecule`) is the
    electron pseudo-species. Mirrors the defensive access of
    :func:`rmgpy.electron_balance.get_electron_species`.
    """
    try:
        return bool(spc.is_electron())
    except (AttributeError, IndexError):
        return False


def _declares_net_electron_count(kinetics):
    """
    ``True`` if ``kinetics`` is one of the plasma rate laws whose ``electrons``
    field carries the same signed net convention as ``Reaction.electrons``.

    Matched on class name against :data:`_NET_ELECTRON_KINETICS_CLASSES` rather
    than by import, because the kinetics package imports this module's
    neighbours; a name comparison keeps the list closed without adding a cycle.
    Subclasses are matched too, so a specialisation of a declared rate law does
    not silently fall back to the blanket refusal.
    """
    for cls in type(kinetics).__mro__:
        if cls.__name__ in _NET_ELECTRON_KINETICS_CLASSES:
            return True
    return False


def _place_declared_electrons(reaction, electron, reactant_count, product_count):
    """
    The list-append primitive: return copies of ``reaction.reactants`` and
    ``reaction.products`` with ``reactant_count`` and ``product_count`` copies of
    ``electron`` appended to the respective side. ``reaction`` is not mutated.

    This is deliberately a SECOND primitive rather than a widening of
    :func:`rmgpy.electron_balance.expand_electrons`. That function is net-derived
    by contract — it appends ``-electrons`` electrons to the reactants and
    nothing to the products — and it is shared with the export path
    (``chemkin.pyx``, both Cantera YAML writers, ``Reaction.is_balanced``), none
    of which has a family declaration in hand to widen it with. Widening it would
    have meant either threading a declaration through four export call sites or
    giving it a default that silently reproduced today's behaviour, and the
    export writers' whole purpose is to be the one place the net scalar IS
    authoritative.

    The two stay in step by construction and by test: on every shape
    ``expand_electrons`` can express — all the one-sided ones — this function
    returns exactly the same lists, which
    ``test/rmgpy/i113IonisationPlacementTest.py`` pins for consumption and
    production alike. They also stay in step downstream, because the view this
    one builds is verified with ``get_species_electron_count`` (step 10), the
    very function the export boundary balances with.

    Unlike ``expand_electrons`` this primitive never looks at
    ``reaction.electrons``: the counts come from the caller, which got them from
    the family declaration. That is what keeps the net scalar out of the order
    path.
    """
    reactants = list(reaction.reactants)
    products = list(reaction.products)
    reactants.extend([electron] * reactant_count)
    products.extend([electron] * product_count)
    return reactants, products


def resolve_electron_placement(reaction, species_list):
    """
    Resolve a family-declared, reactor-facing electron-placement view of
    `reaction`.

    `reaction` is a canonical (database-representation) reaction whose electron
    is carried as ``electrons`` metadata; `species_list` is the model's species
    list, which must contain exactly one canonical electron species.

    Returns a new :class:`rmgpy.reaction.Reaction` — the view — with copied
    participant lists, the family-declared number of electrons inserted on each
    side (one reactant electron and no product electron, for attachment; one and
    two, for an ionisation family), and ``electrons = 0``. The view shares the
    species and kinetics *objects* with the canonical reaction by reference
    (the reactor indexes species by identity), but the canonical reaction
    itself is never mutated: same ``electrons``, same participant lists, after
    this call as before. The view is produced in the reaction's own stored
    orientation, and is a reactor-facing object, never to be written back into
    any database structure.

    Raises :class:`~rmgpy.exceptions.ElectronPlacementError`, always naming the
    offending family/reaction, unless ALL of the following hold: the family
    carries a well-formed placement declaration in
    :data:`FAMILY_ELECTRON_PLACEMENT`; ``reaction.electrons`` equals the declared
    net change, ``product_count - reactant_count``; no explicit electron already
    appears among the participants; the kinetics are a currently supported,
    non-pressure-dependent form whose own electron count, if it has one at all,
    is the declared net change; exactly one canonical electron species is
    resolvable from `species_list`; the finished view balances in the ``E``
    pseudo-element, which is what proves the declared sides were the correct
    sides; and the plasma reaction order implied by the rate coefficient's
    dimensionality, where it is resolvable at all, equals the number of reactants
    in the finished view.

    That last one is a CROSS-CHECK, never an order source: placement is not
    derived from kinetics units here or anywhere. ``get_plasma_rate_order``
    returns ``None`` for any rate law outside its units table, and that silence
    is recorded as :data:`RATE_ORDER_UNRESOLVABLE` — its own outcome — rather
    than being allowed to read as agreement.

    ``is_forward`` and ``reversible`` are NOT acceptance conditions — see the
    module docstring and step 6 for why, and for where the reversibility
    protection lives instead.
    """
    # 1. Family attribution. Placement is family-declared; a reaction that
    #    carries no family cannot name a declaration and must not fall back to
    #    its net electron count.
    family = getattr(reaction, 'family', None)
    if not family:
        raise ElectronPlacementError(
            'Reaction {0!s} carries no family attribution; electron placement is '
            'family-declared and cannot be inferred from the net electron count '
            '(Reaction.electrons).'.format(reaction))

    net_electrons = getattr(reaction, 'electrons', 0) or 0
    n_explicit_reactants = sum(1 for spc in (reaction.reactants or []) if _is_electron(spc))
    n_explicit_products = sum(1 for spc in (reaction.products or []) if _is_electron(spc))
    n_explicit = n_explicit_reactants + n_explicit_products

    # 2. Double representation is fatal, whatever the family: an explicit
    #    electron participant AND a nonzero metadata count would double-count
    #    rate order or electron stoichiometry, and silently preferring either
    #    source would hide the corruption.
    if n_explicit and net_electrons != 0:
        raise ElectronPlacementError(
            'Reaction {0!s} (family {1!r}) represents the electron twice: {2:d} explicit '
            'electron participant(s) AND a nonzero metadata electron count '
            '(electrons={3:d}). This double representation would double-count electron '
            'stoichiometry or rate order; refusing to prefer either source.'.format(
                reaction, family, n_explicit, net_electrons))

    # 3. Family placement declaration. Absence is a named failure, never a
    #    fallback to net-derived inference.
    declaration = FAMILY_ELECTRON_PLACEMENT.get(family)
    if declaration is None:
        raise ElectronPlacementError(
            'Family {0!r} has no electron-placement declaration (reaction {1!s}, '
            'electrons={2:d}); refusing to infer electron placement from the net '
            'electron count. Only families declared in FAMILY_ELECTRON_PLACEMENT '
            'can resolve.'.format(family, reaction, net_electrons))
    #    The shape check is a general one: any pair of non-negative integers
    #    that places at least one electron somewhere is a declaration this
    #    resolver can honour. What it refuses is a MALFORMED declaration —
    #    including the superseded ``(side, count)`` spelling, which would
    #    otherwise be unpacked as a count pair and place electrons by accident.
    if (not isinstance(declaration, tuple) or len(declaration) != 2
            or not all(isinstance(n, int) and not isinstance(n, bool) and n >= 0
                       for n in declaration)
            or declaration == (0, 0)):
        raise ElectronPlacementError(
            'Family {0!r} carries the malformed electron-placement declaration '
            '{1!r} (reaction {2!s}). A declaration must be a '
            '(reactant_count, product_count) pair of non-negative integers '
            'placing at least one electron; refusing to interpret it.'.format(
                family, declaration, reaction))
    reactant_count, product_count = declaration

    # 4. A reaction that already carries its electron explicitly (with a zero
    #    metadata count) is already in reactor form; there is nothing for this
    #    resolver to place, and re-placing would add a second electron.
    if n_explicit:
        raise ElectronPlacementError(
            'Reaction {0!s} (family {1!r}) already carries {2:d} explicit electron '
            'participant(s); the placement view would add another. Nothing to '
            'resolve.'.format(reaction, family, n_explicit))

    # 5. Validate the metadata count AGAINST the declaration. The count is
    #    never the order source; a count that contradicts the declaration is a
    #    reaction of a different shape handed to the wrong resolver.
    #
    #    This generalises from ``electrons == -count`` with NO change of
    #    character. The declaration owns both numbers; the scalar is compared
    #    against the one of them it is able to speak to — the net change — and is
    #    read for nothing else. In particular the scalar never says how many
    #    electrons are incident, which is the whole reason the declaration needs
    #    two sides: ``electrons = +1`` is equally consistent with
    #    ``Li + e- -> Li+ + 2e-`` (order 2) and ``Li -> Li+ + e-`` (order 1), and
    #    only the declaration distinguishes them.
    expected_net = product_count - reactant_count
    if net_electrons != expected_net:
        if net_electrons > 0:
            shape = 'ionization-shaped (net electron production)'
        elif net_electrons == 0:
            shape = 'excitation-shaped (no net electron change)'
        else:
            shape = 'attachment-shaped (net electron consumption)'
        raise ElectronPlacementError(
            'Reaction {0!s} carries electrons={1:d}, which is {2}; family {3!r} '
            'declares net {4:d} ({5:d} electron(s) on the reactant side, {6:d} on '
            'the product side). No placement view is defined for this '
            'shape.'.format(reaction, net_electrons, shape, family, expected_net,
                            reactant_count, product_count))

    # 6. Directionality is NOT read from ``is_forward``, and reversibility is
    #    not this module's business. Both were checked here until I-088; the
    #    reasoning that removed them, and where each protection now lives:
    #
    #    ``is_forward`` — the old check refused ``is_forward=False`` on the
    #    premise that a reverse-generated reaction "would put the electron on the
    #    wrong side". It does not. ``KineticsFamily._create_reaction``
    #    (family.py, see its docstring) stores the reaction in FAMILY-FORWARD
    #    molecular orientation in BOTH generation directions — the reverse branch
    #    swaps the reverse-matched lists back — and sets ``is_forward`` only to
    #    record how the match was found. The engine therefore never hands this
    #    resolver a reverse-ORIENTED reaction: for every one of them the
    #    reactant side is the family reactant side, which is exactly the side the
    #    declaration names. The genuinely reverse-oriented objects are built
    #    elsewhere (the training-set flip in ``get_training_set(get_reverse=True)``,
    #    family.py, which sets ``is_forward = False`` AND negates ``electrons``);
    #    those are plain ``Reaction`` objects with no family attribution, refused
    #    at step 1, and their positive ``electrons`` is refused at step 5 besides.
    #    So ``is_forward`` is uninformative about placement, and refusing on it
    #    excluded precisely the chemistry the plasma decks generate.
    #
    #    What the old check was proxying for — "is the electron really a reactant
    #    in the stored orientation?" — is now VERIFIED rather than trusted, at
    #    step 10, from the participants themselves. That check is strictly
    #    stronger: it does not depend on ``is_forward``, on the producer honouring
    #    the orientation invariant, or on ``electrons`` being signed correctly.
    #
    #    ``reversible`` — the old check refused ``reversible=True`` because "a
    #    reversible view would leave the reverse rate of an electron-containing
    #    reaction implicitly defined". That hazard is real, but it is a REVERSE-RATE
    #    policy, not a placement question: a reversible reaction still has a
    #    definite reactant side in its stored orientation, so which side the
    #    electron goes on is not ambiguous. The policy is owned and enforced by the
    #    consumer, ``PlasmaReactor``, which refuses a reversible electron-containing
    #    reaction by name in ``_validate_reactions`` (kr = kf/Keq(Tgas) would price
    #    the electron's thermochemistry at the gas temperature) and again
    #    defensively in ``generate_rate_coefficients``. Both run on the VIEW, after
    #    placement, and the view carries ``reversible`` through unchanged below, so
    #    placement cannot launder it. Keeping a copy here refused the reaction one
    #    step earlier under a misattributed message: the declaration is
    #    ``('reactants', 1)`` and says nothing about reversibility.

    # 7. Kinetics form. The view is reactor-facing, so kinetics the reactor
    #    path cannot support are refused here, by name, rather than surfacing
    #    downstream as an unnamed failure.
    kinetics = getattr(reaction, 'kinetics', None)
    if kinetics is None:
        raise ElectronPlacementError(
            'Reaction {0!s} (family {1!r}) has no kinetics; a placement view '
            'cannot be certified for an unknown rate form.'.format(reaction, family))
    if kinetics.is_pressure_dependent():
        raise ElectronPlacementError(
            'Reaction {0!s} (family {1!r}) carries pressure-dependent kinetics '
            '{2}; pressure-dependent forms are not a supported shape for a '
            'placement view.'.format(
                reaction, family, kinetics.__class__.__name__))
    #
    #    A kinetics-level electron count is normally a SECOND placement source
    #    and fatal. The two plasma rate laws are the exception, and only because
    #    they say so themselves: ``BadnellRRArrhenius.electrons`` and
    #    ``VoronovEIArrhenius.electrons`` are documented as "the net change in
    #    free electron number, signed to the reaction as written", i.e. the same
    #    quantity as ``Reaction.electrons``. For those, the count is validated
    #    against the declaration exactly as the reaction scalar was at step 5 —
    #    a third source that must AGREE, not a second one that must be absent.
    #    Everything else (the charge-transfer forms, whose ``electrons`` counts
    #    electrons transferred, a different quantity no declaration here speaks
    #    to) stays refused, unchanged.
    kinetics_electrons = getattr(kinetics, 'electrons', None)
    if kinetics_electrons and not _declares_net_electron_count(kinetics):
        raise ElectronPlacementError(
            'Reaction {0!s} (family {1!r}) carries kinetics {2} with its own '
            'electron count (electrons={3!r}); a kinetics-level electron field is a '
            'second placement source and would double-represent the '
            'electron.'.format(reaction, family, kinetics.__class__.__name__,
                               kinetics_electrons))
    if _declares_net_electron_count(kinetics):
        try:
            declared_by_kinetics = int(round(kinetics.electrons.value_si))
        except AttributeError:
            declared_by_kinetics = int(kinetics.electrons)
        if declared_by_kinetics != expected_net:
            raise ElectronPlacementError(
                'Reaction {0!s} (family {1!r}) carries kinetics {2} declaring a net '
                'electron change of {3:d}, which disagrees with the family '
                'declaration {4!r} (net {5:d}). The rate law and the family must '
                'describe the same chemistry; refusing to prefer either '
                'source.'.format(reaction, family, kinetics.__class__.__name__,
                                 declared_by_kinetics, declaration, expected_net))

    # 8. Exactly one canonical electron species must be resolvable.
    electrons_in_model = [spc for spc in (species_list or []) if _is_electron(spc)]
    if len(electrons_in_model) == 0:
        raise ElectronPlacementError(
            'No electron species is resolvable from the supplied species list; '
            'reaction {0!s} (family {1!r}) cannot be given an explicit electron '
            'reactant. Add the canonical electron species to the model.'.format(
                reaction, family))
    if len(electrons_in_model) > 1:
        raise ElectronPlacementError(
            'Found {0:d} electron species in the supplied species list; the '
            'canonical electron must be unique for reaction {1!s} (family {2!r}) '
            'to resolve unambiguously.'.format(
                len(electrons_in_model), reaction, family))
    electron = electrons_in_model[0]

    # 9. Build the view. The placement is now authorised by the family
    #    declaration; _place_declared_electrons is used purely as the list-append
    #    primitive (it copies both lists and appends the declared counts, and
    #    never reads the net scalar), and its output is checked against the
    #    declaration. See that function for why it is a second primitive rather
    #    than a widening of the export path's ``expand_electrons``.
    view_reactants, view_products = _place_declared_electrons(
        reaction, electron, reactant_count, product_count)

    placed_reactants = sum(1 for spc in view_reactants if _is_electron(spc))
    placed_products = sum(1 for spc in view_products if _is_electron(spc))
    if placed_reactants != reactant_count or placed_products != product_count:
        # Defensive: the primitive disagreeing with the validated declaration
        # would mean an electron arrived from somewhere other than the
        # declaration, which cannot be repaired here.
        raise ElectronPlacementError(
            'Electron placement for reaction {0!s} (family {1!r}) produced '
            '{2:d} reactant-side and {3:d} product-side electrons, but the family '
            'declares {4:d} and {5:d}; refusing the inconsistent '
            'view.'.format(reaction, family, placed_reactants, placed_products,
                           reactant_count, product_count))

    # 10. Structural verification that the declared side is the RIGHT side.
    #     Steps 3-9 establish only that the placement matches the family's
    #     declaration and the reaction's metadata count. Both of those are
    #     assertions by a producer; this step checks the result against the
    #     participants themselves, which no producer controls.
    #
    #     With n electrons on one side, the ``E`` pseudo-element balances for
    #     exactly one of the two sides: moving them across changes the imbalance
    #     by 2n, which is nonzero for n >= 1. So an E-balanced view is proof that
    #     the electron went to the correct side, independent of ``is_forward``, of
    #     ``electrons``, and of the orientation invariant the generating family
    #     is supposed to honour. This is what replaced the ``is_forward`` refusal
    #     in step 6 -- verification of the outcome in place of a proxy for it, and
    #     it is the check that would have caught the I-086 sign inversion at this
    #     boundary rather than at the Chemkin writer.
    #
    #     Counted with the writers' own rule (``get_species_electron_count``,
    #     imported read-only) so the reactor boundary and the export boundary
    #     cannot disagree about what balances.
    try:
        reactant_e = sum(get_species_electron_count(spc) for spc in view_reactants)
        product_e = sum(get_species_electron_count(spc) for spc in view_products)
    except (AttributeError, IndexError) as exc:
        raise ElectronPlacementError(
            'Electron placement for reaction {0!s} (family {1!r}) cannot be '
            'verified: the E pseudo-element count is undefined for at least one '
            'participant ({2!s}). A view that cannot be checked is refused rather '
            'than trusted.'.format(reaction, family, exc))
    if reactant_e != product_e:
        raise ElectronPlacementError(
            'Electron placement for reaction {0!s} (family {1!r}) put {2:d} '
            'electron(s) on the reactant side and {3:d} on the product side per the '
            'family declaration, but the resulting view does not balance in the E '
            'pseudo-element (E={4:d} on the left, E={5:d} on the right). The '
            'declared sides cannot be the correct sides for this reaction; refusing '
            'the view rather than handing the reactor an unbalanced one.'.format(
                reaction, family, reactant_count, product_count, reactant_e, product_e))

    # 11. Rate-order CROSS-CHECK. A check, never an order source: deriving
    #     placement from kinetics units is what the module docstring forbids, and
    #     this runs on a view the declaration has already fixed. It converts what
    #     would otherwise be a MechanismWriterError at export
    #     (``check_electron_reactant_order``) into a placement-time refusal,
    #     which is strictly earlier and cheaper.
    #
    #     The number to compare against is the number of REACTANTS in the view,
    #     not the number of reactant electrons: the rate coefficient's
    #     dimensionality fixes the overall order, and the electron is one
    #     reactant among them. ``Li + e- -> Li+ + 2e-`` is order 2 with two
    #     reactants; a ``(0, 1)`` declaration for the same net +1 would give one
    #     reactant against the same order-2 coefficient, which is exactly the
    #     factor-of-electron-density error this catches.
    #
    #     ``get_plasma_rate_order`` returns ``None`` for any rate law whose units
    #     are not in its table — a silent gap, not a named failure. A silent
    #     ``None`` must not read as agreement, so it gets its own outcome.
    resolved_order = get_plasma_rate_order(kinetics)
    if resolved_order is None:
        order_outcome = RATE_ORDER_UNRESOLVABLE
        logging.debug(
            'Electron placement for reaction %s (family %r): %s for kinetics %s; '
            'the placement view stands on the family declaration and the E '
            'pseudo-element balance alone.',
            reaction, family, RATE_ORDER_UNRESOLVABLE, kinetics.__class__.__name__)
    elif resolved_order != len(view_reactants):
        raise ElectronPlacementError(
            'Electron placement for reaction {0!s} (family {1!r}) built a view with '
            '{2:d} reactant(s), but its kinetics {3} has a rate coefficient of order '
            '{4:d}. The family declaration {5!r} and the rate law disagree about the '
            'incident-electron order, so the view would be evaluated at the wrong '
            'order — the rate wrong by a factor of the electron density while the '
            'reaction looks well formed. Refusing the view; the declaration is not '
            'corrected from the rate units.'.format(
                reaction, family, len(view_reactants), kinetics.__class__.__name__,
                resolved_order, declaration))
    else:
        order_outcome = '{0} (order {1:d})'.format(RATE_ORDER_AGREES, resolved_order)

    # The view: a plain reactor-facing Reaction. Species and kinetics objects
    # are shared by reference (the reactor indexes species by identity);
    # degeneracy passes through the constructor, which sets the backing field
    # directly and therefore never touches the shared kinetics object.
    view = Reaction(
        label=reaction.label,
        reactants=view_reactants,
        products=view_products,
        kinetics=reaction.kinetics,
        reversible=reaction.reversible,
        duplicate=reaction.duplicate,
        degeneracy=reaction.degeneracy,
        electrons=0,
        is_forward=getattr(reaction, 'is_forward', None),
        # The cross-check's outcome travels with the view, so a reader can tell
        # "the rate order agrees" from "the rate order could not be resolved"
        # without re-deriving either.
        comment='Electron-placement view (family {0!r}, {1}) of: {2!s}'.format(
            family, order_outcome, reaction),
    )
    return view
