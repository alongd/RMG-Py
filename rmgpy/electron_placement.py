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

Scope: non-dissociative electron attachment
(``Plasma_Electron_Attachment``, ``A + e- -> A-``) only. Every other family,
reaction shape, direction, or kinetics form fails by name with
:class:`rmgpy.exceptions.ElectronPlacementError` — there is no silent fallback
to net-derived inference, and no general mechanism to inherit accidentally.
"""

from rmgpy.electron_balance import expand_electrons
from rmgpy.exceptions import ElectronPlacementError
from rmgpy.reaction import Reaction

__all__ = [
    'FAMILY_ELECTRON_PLACEMENT',
    'resolve_electron_placement',
]

#: Family-level electron placement declarations, keyed by family label.
#: Each value is ``(side, count)``: the side of the reaction the family places
#: its electrons on, in the forward direction, and how many it places there.
#: This mapping — not ``Reaction.electrons`` — is the order source for
#: :func:`resolve_electron_placement`. A family absent from this mapping has NO
#: placement declaration and resolves to a named failure, never to a
#: net-derived guess. Only non-dissociative electron attachment is declared in
#: this increment; ionization- and excitation-type families must gain their own
#: declarations (and their own validation) before they can resolve.
FAMILY_ELECTRON_PLACEMENT = {
    'Plasma_Electron_Attachment': ('reactants', 1),
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
    this call as before. The view is directional — it is produced for, and
    validated in, the forward direction only — and is a reactor-facing object,
    never to be written back into any database structure.

    Raises :class:`~rmgpy.exceptions.ElectronPlacementError`, always naming the
    offending family/reaction, unless ALL of the following hold: the family
    carries a placement declaration (this increment: only
    ``Plasma_Electron_Attachment``); ``reaction.electrons`` equals the declared
    net consumption (-1); no explicit electron already appears among the
    participants; the reaction is irreversible and in its forward direction;
    the kinetics are a currently supported, non-pressure-dependent form; and
    exactly one canonical electron species is resolvable from `species_list`.
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
            '(non-dissociative attachment) is implemented.'.format(family, declaration))

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
            'declares single-electron attachment (net {4:d}). No placement view is '
            'defined for this shape in this increment.'.format(
                reaction, net_electrons, shape, family, expected_net))

    # 6. Directionality. The declaration authorises the FORWARD direction of
    #    an irreversible attachment; a reverse-direction or reversible reaction
    #    would put the electron on the wrong side, or on both.
    if getattr(reaction, 'is_forward', None) is False:
        raise ElectronPlacementError(
            'Reaction {0!s} (family {1!r}) was generated in the reverse direction '
            '(is_forward=False); the electron-placement declaration is defined for '
            'the forward direction only, and a direction-agnostic view is not '
            'safe.'.format(reaction, family))
    if getattr(reaction, 'is_forward', None) is not True:
        # Reaction.__init__ (and TemplateReaction) default is_forward=None, so an
        # unknown direction is NOT the same as a forward one. Placing the
        # electron on the family-declared forward (reactant) side of a reaction
        # whose direction was never established would silently manufacture a
        # forward-direction view from ambiguous input. The pre-integration
        # reactor refused every nonzero-metadata reaction outright; a
        # direction-unknown reaction must be refused here too, by name, rather
        # than resolved and accepted.
        raise ElectronPlacementError(
            'Reaction {0!s} (family {1!r}) has an unspecified reaction direction '
            '(is_forward={2!r}); the electron-placement declaration is defined for '
            'the explicit forward direction only, so a direction-unknown reaction '
            'cannot be given a forward-side electron and is refused rather than '
            'silently accepted.'.format(
                reaction, family, getattr(reaction, 'is_forward', None)))
    if getattr(reaction, 'reversible', False):
        raise ElectronPlacementError(
            'Reaction {0!s} (family {1!r}) is reversible; the placement view is '
            'directional and family {1!r} declares irreversible attachment. A '
            'reversible view would leave the reverse rate of an electron-containing '
            'reaction implicitly defined.'.format(reaction, family))

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
            '{2}; pressure-dependent forms are not a supported shape for the '
            'attachment placement view.'.format(
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
