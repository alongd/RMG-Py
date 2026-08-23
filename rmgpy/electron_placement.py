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
canonical electron species inserted on the side the reaction's *family*
declares, and ``electrons = 0`` in the view. The canonical database reaction is
never modified — it keeps its metadata count and its original participant
lists, so the database representation (and the training-set label convention
built on it) is untouched.

The order source is the family-level placement declaration in
:data:`FAMILY_ELECTRON_PLACEMENT`, established fresh each time the resolver
runs. It is NEVER the net scalar: ``electrons`` is only *validated against* the
declaration, and :func:`rmgpy.electron_balance.expand_electrons` — which is
purely net-derived — is used only as a low-level list-append primitive after
the declaration has authorised the placement, with the result checked against
the declaration afterwards.

Scope: the families named in :data:`FAMILY_ELECTRON_PLACEMENT`, in the
single-electron-on-the-reactant-side shape they declare. Every other family,
reaction shape, or kinetics form fails by name with
:class:`rmgpy.exceptions.ElectronPlacementError` — there is no silent fallback
to net-derived inference, and no general mechanism to inherit accidentally.

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

from rmgpy.electron_balance import expand_electrons, get_species_electron_count
from rmgpy.exceptions import ElectronPlacementError
from rmgpy.reaction import Reaction

__all__ = [
    'FAMILY_ELECTRON_PLACEMENT',
    'resolve_electron_placement',
]

#: Family-level electron placement declarations, keyed by family label.
#: Each value is ``(side, count)``: the side of the reaction the family places
#: its electrons on, in its FAMILY-FORWARD orientation, and how many it places
#: there. This mapping — not ``Reaction.electrons`` — is the order source for
#: :func:`resolve_electron_placement`. A family absent from this mapping has NO
#: placement declaration and resolves to a named failure, never to a
#: net-derived guess.
#:
#: Both declared families place a single electron on the reactant side, for
#: different chemistry: ``Plasma_Electron_Attachment`` is non-dissociative
#: attachment (``A + e- -> A-``); ``Cation_R_Recombination`` is cation-radical
#: recombination (``Li+ + R. + e- -> R-Li``), which the plasma decks reach in the
#: family's REVERSE generation direction — RMG has only the neutral ``R-Li``, so
#: it matches the product template and reconstructs ``Li+ + R.``. The declaration
#: is about the family's forward molecular orientation, which is the orientation
#: the engine stores in BOTH generation directions (see step 6 below), so one
#: entry covers both.
#:
#: Ionization- and excitation-type families (net electron production, or none)
#: must gain their own declarations and their own validation before they can
#: resolve; the shape check in step 5 refuses them by name until then.
FAMILY_ELECTRON_PLACEMENT = {
    'Plasma_Electron_Attachment': ('reactants', 1),
    'Cation_R_Recombination': ('reactants', 1),
}


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


def resolve_electron_placement(reaction, species_list):
    """
    Resolve a family-declared, reactor-facing electron-placement view of
    `reaction`.

    `reaction` is a canonical (database-representation) reaction whose electron
    is carried as ``electrons`` metadata; `species_list` is the model's species
    list, which must contain exactly one canonical electron species.

    Returns a new :class:`rmgpy.reaction.Reaction` — the view — with copied
    participant lists, the electron inserted on the family-declared side
    (reactants, for attachment), and ``electrons = 0``. The view shares the
    species and kinetics *objects* with the canonical reaction by reference
    (the reactor indexes species by identity), but the canonical reaction
    itself is never mutated: same ``electrons``, same participant lists, after
    this call as before. The view is produced in the reaction's own stored
    orientation, and is a reactor-facing object, never to be written back into
    any database structure.

    Raises :class:`~rmgpy.exceptions.ElectronPlacementError`, always naming the
    offending family/reaction, unless ALL of the following hold: the family
    carries a placement declaration in :data:`FAMILY_ELECTRON_PLACEMENT`;
    ``reaction.electrons`` equals the declared net consumption (-1); no explicit
    electron already appears among the participants; the kinetics are a
    currently supported, non-pressure-dependent form carrying no electron count
    of their own; exactly one canonical electron species is resolvable from
    `species_list`; and the finished view balances in the ``E`` pseudo-element,
    which is what proves the declared side was the correct side.

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
    side, declared_count = declaration
    if side != 'reactants' or declared_count != 1:
        # Defensive: nothing else is implemented in this increment.
        raise ElectronPlacementError(
            'Family {0!r} declares electron placement {1!r}, which this resolver '
            'does not support; only a single electron on the reactant side '
            'is implemented.'.format(family, declaration))

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
    expected_net = -declared_count
    if net_electrons != expected_net:
        if net_electrons > 0:
            shape = 'ionization-shaped (net electron production)'
        elif net_electrons == 0:
            shape = 'excitation-shaped (no net electron change)'
        else:
            shape = 'multi-electron attachment-shaped'
        raise ElectronPlacementError(
            'Reaction {0!s} carries electrons={1:d}, which is {2}; family {3!r} '
            'declares single-electron consumption (net {4:d}). No placement view is '
            'defined for this shape in this increment.'.format(
                reaction, net_electrons, shape, family, expected_net))

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
    kinetics_electrons = getattr(kinetics, 'electrons', None)
    if kinetics_electrons:
        raise ElectronPlacementError(
            'Reaction {0!s} (family {1!r}) carries kinetics {2} with its own '
            'electron count (electrons={3!r}); a kinetics-level electron field is a '
            'second placement source and would double-represent the '
            'electron.'.format(reaction, family, kinetics.__class__.__name__,
                               kinetics_electrons))

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
    #    declaration; expand_electrons is used purely as the list-append
    #    primitive (it copies both lists and appends -electrons electrons to
    #    the reactants), and its output is checked against the declaration.
    view_reactants, view_products = expand_electrons(reaction, [electron])

    placed_reactants = sum(1 for spc in view_reactants if _is_electron(spc))
    placed_products = sum(1 for spc in view_products if _is_electron(spc))
    if placed_reactants != declared_count or placed_products != 0:
        # Defensive: expand_electrons disagreeing with the validated
        # declaration would mean the net scalar and the declaration diverged
        # after validation, which cannot be repaired here.
        raise ElectronPlacementError(
            'Electron placement for reaction {0!s} (family {1!r}) produced '
            '{2:d} reactant-side and {3:d} product-side electrons, but the family '
            'declares {4:d} on the reactant side; refusing the inconsistent '
            'view.'.format(reaction, family, placed_reactants, placed_products,
                           declared_count))

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
            'electron(s) on the {3} side per the family declaration, but the '
            'resulting view does not balance in the E pseudo-element (E={4:d} on '
            'the left, E={5:d} on the right). The declared side cannot be the '
            'correct side for this reaction; refusing the view rather than '
            'handing the reactor an unbalanced one.'.format(
                reaction, family, declared_count, side, reactant_e, product_e))

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
        comment='Electron-placement view (family {0!r}) of: {1!s}'.format(
            family, reaction),
    )
    return view
