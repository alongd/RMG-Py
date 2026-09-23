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
I-134 -- reaction identity includes the electron, per side.

The defect these tests pin: RMG decided whether a proposed reaction was already
in the model by comparing reactant and product *references* and *labels*, and
nothing else. The canonical representation keeps the electron out of the
participant lists and in the scalar ``Reaction.electrons``, so two reactions over
the same heavy species were indistinguishable to that comparison. The two shipped
lithium plasma channels are exactly that shape --

    electron-impact ionisation   Li  + e-  =>  Li+ + 2 e-    placement (1, 2)
    radiative recombination      Li+ + e-  =>  Li  + hv      placement (1, 0)

-- so the model kept whichever was offered first and discarded the other as a
duplicate, silently, at ``INFO`` level, with the run exiting 0. A deck asking for
both wrote a mechanism with a cation source and no cation loss channel.

**Why the net scalar is not the fix, which is the subtle part.** The two carry
``electrons = +1`` and ``-1``: exactly equal and opposite, which is the relation a
genuine reverse pair has. A comparison built on the net count therefore also calls
them the same reaction -- and ``Reaction.is_isomorphic`` did, measured, before this
change. They are not reverses: reversing the ionisation gives three-body
recombination, ``(2, 1)``, third order, not the radiative channel's ``(1, 0)``,
second order. Only the per-side placement separates them, which is why the repair
runs through :func:`rmgpy.electron_balance.get_electron_placement_counts` and not
through ``reaction.electrons``.

**The blast radius, which is why the repair is in the identity predicates and not
in the plasma code.** The blind comparison is shared by all four return sites of
``check_for_existing_reaction`` -- two of them the family branches -- and by
``Reaction.is_isomorphic``, which ``ReactionModel.merge``,
``KineticsLibrary.check_for_duplicates`` and the pressure-dependence code all
consume. ``TestChargedFamilyReactionsAreNotConfused`` is the proof that this was
never a lithium special case: it is a real family-generated pair, out of the
database, with no lithium in it.

The classes below pin, in order: the contract for the overwhelming majority of
reactions, which have no electron placement declaration and whose verdicts must
not move at all; that genuine duplicates still collapse, which is what the check
is for; that a real charged family's reactions are no longer confused with each
other; that the real lithium mechanism reaches the model with both channels; and
that a kinetics library carrying both channels can be loaded, which it could not
be before.
"""

import hashlib
import math
import os
import re
import shutil

import numpy as np

import pytest

import rmgpy.data.rmg
from rmgpy import constants, settings
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import KineticsLibrary, LibraryReaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.electron_balance import (get_electron_placement_counts, get_placement_declaration,
                                    get_plasma_rate_order)
from rmgpy.electron_placement import FAMILY_ELECTRON_PLACEMENT, resolve_electron_placement
from rmgpy.exceptions import DatabaseError, ElectronPlacementError
from rmgpy.kinetics import (Arrhenius, BadnellRRArrhenius, Chebyshev,
                            ElectronCollisionPlasma, Lindemann, MultiArrhenius,
                            MultiPDepArrhenius, PDepArrhenius, ThirdBody,
                            TwoTemperaturePlasma, VoronovEIArrhenius)
from rmgpy.kinetics.model import get_reaction_order_from_rate_coefficient_units
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import CoreEdgeReactionModel, are_identical_species_references
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.species import Species

IONISATION = 'PlasmaElectronImpactIonization'
RECOMBINATION = 'PlasmaRadiativeRecombination'
ATTACHMENT = 'Plasma_Electron_Attachment'

#: Triplet dioxygen, the reactant the attachment family is trained on.
O2_ADJACENCY_LIST = """
multiplicity 3
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
"""


def _net_rule_counts(electrons):
    """The placement the net-derived rule produces: the rule
    :func:`rmgpy.electron_balance.expand_electrons` falls back to, restated here
    independently so the tests below compare against a second statement of it
    rather than against the implementation under test."""
    if electrons < 0:
        return -electrons, 0
    return 0, electrons


def _species(smiles):
    return Species(molecule=[Molecule().from_smiles(smiles)])


# ---------------------------------------------------------------------------
# Selecting the reaction a check is about, out of a library the check does not
# own.
#
# The checks below once reached the two lithium channels by asserting that the
# library holding them had exactly one entry and taking ``reactions[0]``. That
# couples every one of them to every unrelated addition anywhere in that
# library: when the database grew an argon radiative-recombination entry, four
# checks went red and a fifth was disabled, and none of what they guard had
# moved. Worse, had the argon entry landed *first* in the list, the checks would
# have gone on passing while measuring argon.
#
# The first replacement named the reaction it meant by its PHYSICAL ATTRIBUTES:
# the participant graphs, then -- as each round found a substitution the current
# key accepted -- the rate order, then the sense of the Te dependence. That
# converged on nothing. Three rounds running, the bypass simply moved to an
# attribute the key did not yet read, and it always can: the attribute space is
# open, so for any finite key someone can author a rate law that agrees with the
# real entry on every attribute in it. A selector built that way is chasing a
# receding target, and worse, each attribute it gains is an attribute on which a
# wrong value now SILENTLY REDIRECTS THE SELECTION instead of failing.
#
# So this round inverts it. What a check means by "the lithium radiative
# recombination" is not "a second-order reaction over Li+ and Li whose rate falls
# with Te". It is "the entry ``[Lip] => [Li]`` of the library
# ``PlasmaRadiativeRecombination``" -- an identity, and identity is carried by
# PROVENANCE, which is complete by construction: two entries of one library, or
# two entries of different libraries, are different entries whatever their
# coefficients say, and no rate law anyone can author changes which entry it is
# written in.
#
# The rate-order and Te-response work is not deleted. It moves to the other side
# of the selection, where it belongs: the entry is selected by provenance and
# then ASSERTED to carry the order, the response, the fit and the reversibility
# its channel implies. A wrong value is then a red check naming the attribute
# that is wrong, instead of a quiet change of subject.
# ---------------------------------------------------------------------------

#: Adjacency lists for the plasma species the entries below are ASSERTED to be over.
#: They are not part of the selection key -- provenance is -- so a species that
#: fails to build here is a red assertion and never a redirected selection.
#:
#: Adjacency lists and deliberately NOT SMILES. Measured on this tree:
#: ``Molecule().from_smiles('[Ar+]')`` returns Ar(2+) -- net charge +2, closed
#: shell, no unpaired electron -- and is therefore *not* isomorphic to the
#: ``[Arp]`` the shipped library holds. ``[He+]`` is lost the same way. The
#: radical cations are the ones that break, because the parser reads the radical
#: electron as a second unit of charge; ``[Li+]``, closed-shell, happens to
#: survive. An adjacency list states unpaired electrons and formal charge
#: separately, and round-trips all four species below exactly --
#: ``TestTheSelectorFailsLoudly`` carries the measurement.
LITHIUM_ATOM = '1 Li u1 p0 c0'
LITHIUM_CATION = '1 Li u0 p0 c+1'
ARGON_ATOM = '1 Ar u0 p4 c0'
ARGON_CATION = '1 Ar u1 p3 c+1'


def _structure(adjacency_list, label=''):
    """A one-species :class:`Species` from an adjacency list."""
    return Species(label=label, molecule=[Molecule().from_adjacency_list(adjacency_list)])


def _same_multiset(found, wanted):
    """True when two species lists are the same multiset under graph isomorphism.

    Greedy consumption is exact here, rather than merely convenient, because
    isomorphism is an equivalence relation: if a wanted structure matches any
    unconsumed candidate then it matches every candidate in that same class, so
    no pairing choice can strand a later structure that a different choice would
    have satisfied.
    """
    if len(found) != len(wanted):
        return False
    remaining = list(found)
    for query in wanted:
        for position, candidate in enumerate(remaining):
            if candidate.is_isomorphic(query):
                del remaining[position]
                break
        else:
            return False
    return True


#: The kinds of "no answer" :func:`rate_order` and
#: :func:`electron_temperature_response` can return. They are the ``kind`` field
#: of an :class:`Undetermined`, not prefixes of a message: a caller branches on
#: the kind and never parses prose.
NO_KINETICS = 'no-kinetics'
UNREADABLE_UNITS = 'unreadable-units'
UNSUPPORTED_CLASS = 'unsupported-class'
INCONSISTENT_COMPONENTS = 'inconsistent-components'
EMPTY_COMPOSITE = 'empty-composite'
COMPONENT_WITHOUT_UNITS = 'component-without-units'
PRESSURE_DEPENDENT_ORDER = 'pressure-dependent-order'
NO_TE_WINDOW = 'no-te-window'
NOT_EVALUABLE = 'not-evaluable'
NOT_POSITIVE = 'not-positive'
FLAT_IN_TE = 'flat-in-te'


class Undetermined:
    """Why a rate coefficient states no order, or no direction in Te.

    A **value with a machine-readable ``kind``**, not a diagnostic string. The
    previous design returned prose and distinguished "unsupported class" from
    "unreadable units" from "no kinetics" only by the first word of that prose, so
    a caller that wanted to branch on the kind had to parse English -- and any
    rewording of a message silently changed the behaviour of every such caller.
    The kind is now the thing; the ``detail`` exists only to be read by a human in
    a failure message.

    Equality is on the kind alone, so ``rate_order(k) == UNSUPPORTED_CLASS`` is
    not how this is used -- ``rate_order(k).kind is UNSUPPORTED_CLASS`` is. Two
    undetermined results of the same kind compare equal even when their details
    differ, because the kind is what a caller acts on. Against an ``int`` or
    against :data:`RISES` / :data:`FALLS` it compares unequal, which is what makes
    it safe to hand straight to an ``== expected`` assertion.
    """

    def __init__(self, kind, detail):
        self.kind = kind
        self.detail = detail

    def __eq__(self, other):
        if isinstance(other, Undetermined):
            return self.kind == other.kind
        return NotImplemented

    def __ne__(self, other):
        result = self.__eq__(other)
        return result if result is NotImplemented else not result

    def __hash__(self):
        return hash((Undetermined, self.kind))

    def __str__(self):
        return '{0}: {1}'.format(self.kind, self.detail)

    def __repr__(self):
        return 'Undetermined({0!r}, {1!r})'.format(self.kind, self.detail)


def _identify(species):
    """A species rendered the way production's ``==`` reads it: by object identity.

    Two colliders that print the same and compare unequal are the whole failure
    mode this renders for, so a message that showed only the label would leave the
    reader unable to see why the guard split them.
    """
    return 'None' if species is None else '{0}@{1:#x}'.format(species, id(species))


def cations_produced_but_never_consumed(reactions):
    """The cation species `reactions` makes and no reaction in `reactions` uses up.

    Species are collected and matched by **object identity**, which is how the
    model itself compares them: ``are_identical_species_references`` (``rmgpy/rmg/
    model.py``) compares ``rxn1.reactants == rxn2.reactants``, and a ``Species``
    inherits object identity for ``==``. Until this round the check that consumes
    this reduced each species to ``species.label`` and compared strings, so an
    unrelated reactant that happened to render under a produced cation's label
    satisfied it -- and in this format the label IS a rendering of the chemistry
    (``[Lip] => [Li]``), so two different species sharing one is not a contrivance,
    it is the normal case for a species reconstructed on the other side of a
    boundary. A check that a cation is consumed must be about the cation, not about
    a string that names it.

    Returns the unconsumed cations in the order they were produced, so a failure
    message can name them.
    """
    produced = {}
    for reaction in reactions:
        for product in reaction.products:
            if product.molecule and product.molecule[0].get_net_charge() > 0:
                produced.setdefault(id(product), product)
    consumed = set()
    for reaction in reactions:
        for reactant in reaction.reactants:
            if id(reactant) in produced:
                consumed.add(id(reactant))
    return [species for key, species in produced.items() if key not in consumed]


def rate_order(kinetics):
    """The total order of `kinetics`' rate coefficient.

    Returns the integer order, or an :class:`Undetermined` whose ``kind`` says
    why there is no integer. Keeping the reasons apart is the point -- "this
    candidate never stated a rate at all", "this candidate stated units nothing
    can parse", "this class states no order by any route" and "this class's
    components disagree about the order" are different situations with different
    repairs, and collapsing them into a bare ``None`` is how a caller ends up
    unable to say what went wrong. The kinds this can return are
    :data:`NO_KINETICS`, :data:`UNREADABLE_UNITS`, :data:`UNSUPPORTED_CLASS`,
    :data:`INCONSISTENT_COMPONENTS`, :data:`EMPTY_COMPOSITE`,
    :data:`COMPONENT_WITHOUT_UNITS` and :data:`PRESSURE_DEPENDENT_ORDER`.

    :data:`EMPTY_COMPOSITE` is the round-84 split: a ``MultiArrhenius`` holding no
    components used to be reported as an *unsupported class*, conflating malformed
    data in a class this reader fully supports with a class it has never heard of.
    The two have different repairs -- fix the entry, versus teach this function a
    class -- and a caller that cannot tell them apart cannot make either.

    **Two routes, in a fixed precedence, and why there are two.** The primary
    route is dimensional: read the order off the rate coefficient's units through
    RMG's own :func:`get_reaction_order_from_rate_coefficient_units`, which reads
    any spelling. Reading it off the units rather than off
    ``len(reaction.reactants)`` is the whole point -- the participant lists do not
    carry the electron, so counting reactants cannot see the electron order,
    which is precisely what separates the two channels this function exists to
    tell apart.

    **Where the units live is class-specific, and a single ``A`` attribute does
    not find them.** ``Chebyshev`` states ``kunits``; ``ThirdBody`` states them on
    ``arrheniusLow``; ``PDepArrhenius``, ``MultiArrhenius`` and
    ``MultiPDepArrhenius`` state them on each component. Until this round all five
    were reported ``unsupported``, which was wrong in the loudest possible way --
    it is the outcome that means "no check happened here", and it was being
    returned for classes that state their dimensions perfectly well. Measured on
    database 96f2afa4a, 16 shipped plasma entries are ``ThirdBody`` alone.
    :func:`_rate_coefficient_units` walks those routes and requires every
    component of a composite to agree; a composite whose components disagree is
    :data:`INCONSISTENT_COMPONENTS`, never a silent pick of the first.

    **Agreeing means agreeing about the ORDER, and every component has to state
    one.** Both halves of that were wrong until this round and both were measured
    on the built module, not read off the code. A component that stated no units
    at all was skipped, so ``MultiArrhenius`` holding one ``Arrhenius`` in
    ``m^3/(mol*s)`` and one stating nothing came back as order 2 -- the valid
    component's answer, wearing the composite's name; that is now
    :data:`COMPONENT_WITHOUT_UNITS`. And agreement was tested on the unit STRING,
    so ``['m^3/(mol*s)', 'cm^3/(mol*s)']`` -- two spellings of order 2 -- was
    reported :data:`INCONSISTENT_COMPONENTS`. Components are now normalised to
    their order before being compared, and the spellings are kept in the message.

    A falloff form with *both* limits (``Lindemann``, ``Troe``) is deliberately
    not given an order: its two limits differ in order by one, by construction, so
    the honest answer is :data:`PRESSURE_DEPENDENT_ORDER` and not either limit.

    The dimensional route is **blind on a whole kinetics class**, and this is
    measured, not assumed: :class:`ElectronCollisionPlasma` stores a tabulated
    cross-section and has no ``A`` attribute at all (``hasattr(k, 'A')`` is
    ``False``), so the units route reported ``unsupported`` for it -- while
    production, in ``rmgpy.electron_balance.get_plasma_rate_order``, defines its
    order as exactly 2, because ``<sigma*v>`` is bimolecular by construction.
    Three of the five shipped plasma kinetics libraries (``PlasmaAir``,
    ``PlasmaAlkali``, ``PlasmaArgon``) and two plasma families are written in that
    class, so the blind spot is not a corner.

    So the fallback route asks the engine. The precedence is deliberate and is
    the opposite of convenience: the units are consulted FIRST so that wherever
    both can speak, this function's answer is an independent reading and not an
    echo. A key that simply called the engine's function would agree with it by
    construction and could not notice the engine's ten-entry unit table going
    stale -- which is what ``test_the_rate_order_read_here_agrees_with_the_engine_s_own_reading``
    exists to catch. The engine is consulted only where the units say nothing,
    which is exactly where there is nothing to be independent about.
    """
    if kinetics is None:
        return Undetermined(NO_KINETICS, 'there is no rate object at all')
    units = _rate_coefficient_units(kinetics)
    if isinstance(units, Undetermined):
        return units
    if units is not None:
        try:
            return get_reaction_order_from_rate_coefficient_units(units)
        except ValueError as error:
            return Undetermined(UNREADABLE_UNITS, '{0} states units {1!r}, which RMG '
                                'cannot read ({2})'.format(type(kinetics).__name__,
                                                           units, error))
    from_engine = get_plasma_rate_order(kinetics)
    if from_engine is not None:
        return from_engine
    return Undetermined(UNSUPPORTED_CLASS,
                        '{0} carries no rate-coefficient units by any route this '
                        'function knows, and the engine states no order for it '
                        'either'.format(type(kinetics).__name__))


#: The classes whose whole content is a list of sub-rates. Named by class and not
#: by the presence of an ``arrhenius`` attribute, because the three do not agree
#: on what an empty one looks like: measured, ``PDepArrhenius()`` holds ``[]``
#: while ``MultiArrhenius()`` and ``MultiPDepArrhenius()`` hold ``None``. Keying on
#: the attribute alone therefore sent two of the three down the "class nobody
#: knows" path and only one down the composite path.
COMPOSITE_CLASSES = (MultiArrhenius, MultiPDepArrhenius, PDepArrhenius)


def _components(kinetics):
    """The sub-rates a composite kinetics object is made of, or ``None``.

    ``[]`` and ``None`` mean different things here: ``[]`` is a composite class
    holding nothing, which is malformed data, and ``None`` is a class that is not
    a composite at all.
    """
    parts = getattr(kinetics, 'arrhenius', None)
    if parts is None:
        return [] if isinstance(kinetics, COMPOSITE_CLASSES) else None
    try:
        return list(parts)
    except TypeError:
        return None


def _rate_coefficient_units(kinetics):
    """The units of `kinetics`' rate coefficient, ``None`` if it states none, or
    an :class:`Undetermined` when it states several that disagree.

    Four routes, tried in order, because RMG's kinetics classes do not agree on
    where the dimensions live. See :func:`rate_order` for what each one covers and
    why a bare ``getattr(kinetics, 'A')`` was not enough.
    """
    direct = getattr(getattr(kinetics, 'A', None), 'units', None)
    if direct is not None:
        return direct

    chebyshev = getattr(kinetics, 'kunits', None)
    if chebyshev is not None:
        return chebyshev

    low = getattr(kinetics, 'arrheniusLow', None) or getattr(kinetics, 'arrhenius_low', None)
    if low is not None:
        high = getattr(kinetics, 'arrheniusHigh', None) or \
            getattr(kinetics, 'arrhenius_high', None)
        if high is not None:
            return Undetermined(
                PRESSURE_DEPENDENT_ORDER,
                '{0} states a low-pressure and a high-pressure limit, whose orders '
                'differ by one by construction, so it has no single order'.format(
                    type(kinetics).__name__))
        return _rate_coefficient_units(low)

    parts = _components(kinetics)
    if parts is not None and not parts:
        return Undetermined(
            EMPTY_COMPOSITE,
            '{0} is a composite rate law holding NO components, so it states no '
            'order because it states no rate. That is a malformed instance of a '
            'class this function fully supports, and it is reported apart from '
            '"unsupported class": the repair is to the entry that authored an '
            'empty composite, not to this reader.'.format(type(kinetics).__name__))
    if parts:
        seen = []
        for position, part in enumerate(parts, start=1):
            units = _rate_coefficient_units(part)
            if isinstance(units, Undetermined):
                return units
            if units is None:
                # A component whose units cannot be read is NOT evidence about the
                # composite's order, and until this round it was silently dropped:
                # the loop skipped it and the composite inherited whatever the
                # readable components said. Measured on the built module, a
                # MultiArrhenius holding one Arrhenius in 'm^3/(mol*s)' and one
                # component stating no units at all came back as order 2 --
                # indistinguishable from the valid component alone. So a shipped
                # entry half of whose rate law is malformed passed the sweep that
                # exists to refuse exactly that, and every order-keyed check
                # downstream read an order that half the rate never stated.
                return Undetermined(
                    COMPONENT_WITHOUT_UNITS,
                    'component {0} of {1} ({2}) states no rate-coefficient units by '
                    'any route this function knows, so the composite has no order '
                    'that is a property of ALL of it. The order the readable '
                    'components state is not the composite\'s order: the repair is '
                    'to the component that states none.'.format(
                        position, type(kinetics).__name__, type(part).__name__))
            if units not in seen:
                seen.append(units)
        if len(seen) > 1:
            # Compare the components by ORDER, not by spelling. Until this round
            # this compared unit STRINGS, so two components of the same dimension
            # written differently were reported as disagreeing: measured,
            # MultiArrhenius['m^3/(mol*s)', 'cm^3/(mol*s)'] came back
            # inconsistent-components while RMG reads both spellings as order 2.
            # "Inconsistent" has to mean the orders disagree, or the outcome that
            # means "no check happened here" fires on a perfectly determinate rate.
            # The spellings stay in the message: a reader of a red run still needs
            # to see what was actually written.
            orders = {}
            for spelling in seen:
                try:
                    order = get_reaction_order_from_rate_coefficient_units(spelling)
                except ValueError as error:
                    return Undetermined(
                        UNREADABLE_UNITS,
                        'component of {0} states units {1!r}, which RMG cannot read '
                        '({2})'.format(type(kinetics).__name__, spelling, error))
                orders.setdefault(order, []).append(spelling)
            if len(orders) > 1:
                return Undetermined(
                    INCONSISTENT_COMPONENTS,
                    '{0} is built from components whose units read as {1} different '
                    'orders ({2}), so it has no single order'.format(
                        type(kinetics).__name__, len(orders),
                        '; '.join('order {0} from {1}'.format(
                            order, ', '.join(repr(u) for u in spellings))
                            for order, spellings in sorted(orders.items()))))
        # `parts` is non-empty and every component either returned above or added
        # its spelling, so `seen` cannot be empty here.
        return seen[0]
    return None


def _describe_order(order):
    return 'rate order {0}'.format(order) if isinstance(order, int) else str(order)


#: The two decisive values of :func:`electron_temperature_response`.
FALLS = 'k falls with Te'
RISES = 'k rises with Te'

#: The electron-temperature band the response is read over, in K. Roughly 0.86 eV
#: to 86 eV: the band a low-temperature plasma model actually runs in, wide enough
#: that a real Te dependence cannot hide inside it and narrow enough to sit inside
#: the declared validity range of every shipped plasma rate but one (Voronov's
#: Tmin is 11 604.5 K, so its window is clamped -- see below).
TE_BAND = (1.0e4, 1.0e6)

#: How far the endpoint ratio must be from 1 before the response is called
#: decisive. Measured margins on the shipped entries, over their own windows:
#: Badnell Li+ 0.0228, Shull-Van-Steenberg Ar+ 0.0499, Voronov Li 390. The nearest
#: is 10x clear of this bar, so the bar is not fitted to the data it passes; and a
#: rate with no Te dependence at all gives exactly 1.0 and is correctly called
#: indecisive rather than being pushed to one side by float noise.
TE_DECISIVE = 2.0

#: The GAS temperature every rate in this file is evaluated at. 298.15 K is the
#: 5 torr argon deck's own ``Tgas``, named in the argon entry's ``longDesc``; it
#: is held fixed while ``Te`` is varied, which is what makes the band below a
#: reading of the *electron*-temperature response and not of a mixture of the two.
#:
#: It is a real parameter of the evaluation and not a formality. The one shipped
#: rate law that is a function of both temperatures, ``TwoTemperaturePlasma``,
#: carries ``exp(Ea_e (Te - T) / (R T Te))``, which is identically 1 when T == Te
#: -- see :func:`evaluate_as_the_solver_does`.
SOLVER_TGAS = 298.15


def evaluate_as_the_solver_does(kinetics, Te, Tgas=SOLVER_TGAS):
    """``k`` as the reactor that runs these kinetics would compute it, at (Tgas, Te).

    This calls ``PlasmaReactor.evaluate_two_temperature_rate_coefficient``
    itself. It is not a restatement of production's dispatch rule, it IS
    production's dispatch rule, reached through a reactor built here for the
    purpose; a second copy of a dispatch table is a second thing to drift.

    WHAT WAS WRONG BEFORE, MEASURED
    -------------------------------
    Until this round the Te-response reading and the magnitude anchors both
    called ``kinetics.get_rate_coefficient(T)``. For three of the four shipped
    plasma rate laws that is the right call and for the fourth it is not, and
    the difference is the whole point of a two-temperature reactor:

    ========================  =========================  ========================
    class                     is a function of           ``get_rate_coefficient(T)``
    ========================  =========================  ========================
    ``VoronovEIArrhenius``    ``Te`` only                is ``k(Te=T)`` -- correct
    ``BadnellRRArrhenius``    ``Te`` only                is ``k(Te=T)`` -- correct
    ``ElectronCollisionPlasma``  ``Te`` only             is ``k(Te=T)`` -- correct
    ``TwoTemperaturePlasma``  ``Tgas`` AND ``Te``        is ``k(T, Te=T)`` -- WRONG
    ========================  =========================  ========================

    The last row collapses the two temperatures onto one, and the Kossyi form's
    electron-activation term ``exp(Ea_e (Te - T) / (R T Te))`` is then exactly 1
    whatever ``Ea_e`` is. Measured on the built module, with the argon entry's own
    ``A``, ``n`` and ``T0``:

        Ea_e = 0 kJ/mol    get_rate_coefficient(20000)               = 144584.52071909094
        Ea_e = 100 kJ/mol  get_rate_coefficient(20000)               = 144584.52071909094
        Ea_e = 0 kJ/mol    two_temp(298.15, 20000)                   = 1.44585e5
        Ea_e = 100 kJ/mol  two_temp(298.15, 20000)                   = 2.61925e22

    Identical to the last digit on the interface the checks used; a factor of
    1.8116e17 apart at the low anchor, and 3.1125e17 at the high one, on the
    interface the solver uses. So ``Ea_e`` was a field of
    the shipped argon entry that could take any value at all with every check in
    this file staying green -- the same defect as the unasserted ``A``, one field
    over. Pinned by ``test_the_one_temperature_interface_cancels_the_electron_
    activation_energy``.

    WHY THE DISPATCH IS NOT A CLASS-NAME TEST
    -----------------------------------------
    ``uses_electron_temperature`` plus the evaluator the class *declares* is what
    ``rmgpy/solver/plasma.pyx`` branches on, and a check keyed on
    ``isinstance(kinetics, TwoTemperaturePlasma)`` would go quietly blind on the
    next Te-dependent class the database learns to author. Delegating to the
    reactor means there is nothing here to keep in step.

    ``Tgas`` is a real argument and not a stand-in: a rate law that is a function
    of both temperatures is being asked what it does at ONE gas temperature, and a
    caller that cares about another one passes it.

    Raises whatever the reactor raises. In particular a kinetics object that does
    NOT declare ``uses_electron_temperature`` is not a function of ``Te`` at all --
    production evaluates it at ``(Tgas, P)`` through the ordinary thermal branch --
    so asking this question of one is a category error, and it surfaces as a
    ``TypeError`` here rather than as a number that looks like an answer.
    """
    if not getattr(kinetics, 'uses_electron_temperature', False):
        raise TypeError(
            '{0} does not declare uses_electron_temperature, so the reactor '
            'evaluates it at the GAS temperature through the ordinary thermal '
            'branch and it is not a function of Te at all. Evaluating it "at Te" '
            'would be reading a gas-temperature rate law off an electron-'
            'temperature axis.'.format(type(kinetics).__name__))
    reactor = PlasmaReactor(T=(Tgas, 'K'), P=(5.0, 'torr'), Te=(Te, 'K'),
                            initial_mole_fractions={})
    return reactor.evaluate_two_temperature_rate_coefficient(kinetics)


def electron_temperature_response(kinetics):
    """Whether `kinetics` rises or falls with electron temperature.

    Returns :data:`RISES`, :data:`FALLS`, or an :class:`Undetermined` whose
    ``kind`` says why the sense could not be read -- :data:`NO_KINETICS`,
    :data:`NO_TE_WINDOW`, :data:`NOT_EVALUABLE`, :data:`NOT_POSITIVE` or
    :data:`FLAT_IN_TE`. A caller branches on the kind and never on the prose.

    **The attribute this reads, and why it is not any of the others.** Two
    reactions that this suite has to tell apart can agree on every other property
    it keys on. Measured, on the shipped database:

    ======================  =========================  =========================
    attribute               lithium radiative recomb.  lithium ionisation rate
                            ``BadnellRRArrhenius``     ``VoronovEIArrhenius``
                            ``(Z=3, N=2)``             ``(Z=3, N=3)``
    ======================  =========================  =========================
    reactant structures     ``[Li+]``                  ``[Li+]``   (same entry)
    product structures      ``[Li]``                   ``[Li]``    (same entry)
    net ``electrons``       ``-1``                     ``-1``      (settable)
    placement counts        ``(1, 0)``                 ``(1, 0)``  (owner-keyed)
    rate order              ``2``                      ``2``
    ``A`` units             ``cm^3/(molecule*s)``      ``cm^3/(molecule*s)``
    **sense of k(Te)**      **falls, ratio 0.0228**    **rises, ratio 390**
    ======================  =========================  =========================

    So pasting the shipped Voronov ionisation fit onto the shipped recombination
    entry -- with ``electrons = -1`` so even the bookkeeping is untouched -- is a
    substitution that changes the reaction and that *every* key in the right-hand
    column above accepts. It is also the mutation closest to the accepted case
    that exists: the replacement rate is not invented, it is another fit shipped
    by this same database for this same element, loaded from the same
    ``voronov.yaml`` table.

    The sense of the Te dependence refuses it, and refuses it for a physical
    reason rather than a bookkeeping one. Electron-impact ionisation is a
    threshold process: below the ionisation energy essentially no electron in the
    distribution can ionise, so k climbs steeply as Te rises through the
    threshold. Radiative recombination is the opposite -- a faster electron is
    less likely to be captured, so k falls monotonically with Te. A rate law that
    falls with Te is not an ionisation rate, whatever it is labelled.

    **It is orthogonal to the rate order, in both directions, and that is
    measured.** The round-76 swap (a third-order ``TwoTemperaturePlasma`` with
    ``n = -4.5`` standing in for three-body recombination) *falls* with Te, ratio
    1e-9 -- so the response alone does not catch it and the order does. The
    Voronov substitution here is second order -- so the order alone does not catch
    it and the response does. Neither key subsumes the other; both are needed, and
    a check that had only ever seen one of the two mutations would have concluded
    the wrong one was sufficient.

    **Why a ratio and not monotonicity.** Measured: the Voronov fit is NOT
    monotonic over its own declared range. ``k`` peaks near 3e5 K and falls again
    above it, because the Voronov form carries a ``U^K exp(-U)`` factor with
    ``U = dE/Te`` that drives ``k`` back to zero as ``Te -> inf``. A monotonicity
    test would therefore have called the shipped ionisation rate unreadable. The
    endpoint ratio over a bounded band is the statement that survives contact with
    the actual functions.

    The window is :data:`TE_BAND` intersected with the rate's own declared
    ``[Tmin, Tmax]``, so a fit is never evaluated outside the range its authors
    stated. The narrowing is real: Voronov's ``Tmin`` of 11 604.5 K raises the low
    end of its window, which is why its ratio is 390 here and 1.3e14 over a band
    that started below its threshold. Either is decisive; the clamped one is the
    one that is honest about the fit.

    **That clamp is also why a declared window has to be asserted separately, and
    it is.** This function reads ``Tmin``/``Tmax`` off the rate law and then
    evaluates inside whatever they say, so narrowing them moves the two points this
    reading is taken at and the reading follows -- a narrowed window can never show
    up here as a wrong response. :meth:`ShippedEntry.assert_as_shipped` therefore
    names the declared window as its own assertion, exactly as it does the pressure
    window. A check that silently adapts to the thing that changed reports nothing.

    **The rate is evaluated through** :func:`evaluate_as_the_solver_does` **, which
    is the reactor's own dispatch.** ``Tgas`` is held at :data:`SOLVER_TGAS` while
    ``Te`` is swept, so what is read is the response to the *electron* temperature.
    Reading it off ``get_rate_coefficient(Te)`` instead -- which is what this did
    until this round -- sweeps both temperatures together on
    ``TwoTemperaturePlasma`` and cancels its ``Ea_e`` term identically.
    """
    if kinetics is None:
        return Undetermined(NO_KINETICS, 'there is no rate object at all')
    low, high = TE_BAND
    declared_low = getattr(getattr(kinetics, 'Tmin', None), 'value_si', None)
    declared_high = getattr(getattr(kinetics, 'Tmax', None), 'value_si', None)
    if declared_low is not None:
        low = max(low, declared_low)
    if declared_high is not None:
        high = min(high, declared_high)
    if not high > low:
        return Undetermined(NO_TE_WINDOW,
                            'no usable Te window on {0}: [{1:g}, {2:g}] K does not '
                            'overlap its declared [{3}, {4}] K'.format(
                                type(kinetics).__name__, TE_BAND[0], TE_BAND[1],
                                declared_low, declared_high))
    try:
        at_low = evaluate_as_the_solver_does(kinetics, low)
        at_high = evaluate_as_the_solver_does(kinetics, high)
    except Exception as error:                                # noqa: BLE001
        return Undetermined(NOT_EVALUABLE,
                            'k(Te) could not be evaluated on {0}: {1}: {2}'.format(
                                type(kinetics).__name__, type(error).__name__, error))
    if at_low <= 0.0 or at_high <= 0.0:
        return Undetermined(NOT_POSITIVE,
                            'k(Te) is not positive on {0}: k({1:g} K) = {2:g}, '
                            'k({3:g} K) = {4:g}'.format(type(kinetics).__name__,
                                                        low, at_low, high, at_high))
    ratio = at_high / at_low
    if ratio >= TE_DECISIVE:
        return RISES
    if ratio <= 1.0 / TE_DECISIVE:
        return FALLS
    return Undetermined(FLAT_IN_TE,
                        'k is flat in Te on {0}: k({1:g} K)/k({2:g} K) = {3:g}, within '
                        'a factor of {4:g} of 1'.format(type(kinetics).__name__, high,
                                                        low, ratio, TE_DECISIVE))


def _describe_response(response):
    return response if isinstance(response, str) else str(response)


# ---------------------------------------------------------------------------
# Rate magnitude, anchored to the publication rather than to the database
# ---------------------------------------------------------------------------
#
# ``shipped_fit`` compares the entry's rate law against another object built by
# the same ``(Z, N)`` lookup into the same shipped YAML table. That catches a fit
# for the wrong element, which is what it was written for, and it cannot catch
# anything the two sides share: corrupt ``badnell.yaml`` and BOTH sides move
# together, and an entry that hand-authors the shipped numbers inline instead of
# looking them up is indistinguishable from one that looks them up. It is also
# silent about magnitude wherever no shipped table exists -- which is the argon
# entry, whose ``A`` could be replaced with ``1e100`` and every check stayed
# green, measured.
#
# The anchors below close both. Each one states a *published closed form* and its
# *published parameters*, transcribed into this file from the paper the entry's
# own ``shortDesc`` cites, and evaluates it here. Nothing on the reference side
# is read from the database, so a corrupted table, a hand-authored inline copy
# and a silently re-scaled ``A`` all move one side only.
#
# WHAT THIS COSTS, ON THE DAY THE DATABASE RE-FITS AN ENTRY LEGITIMATELY
# ----------------------------------------------------------------------
# It costs a red check and a deliberate edit here, naming the new source. That is
# the intended price and not a defect: a re-fit of argon's recombination rate is a
# scientific decision, and the whole complaint this answers is that such a change
# was invisible. The anchor makes it visible exactly once, at the moment it
# happens, to the person making it.
#
# The alternative that was considered and rejected: asserting that the entry was
# built by the shipped-table route rather than from inline numbers. It is not
# available -- measured, ``BadnellRRArrhenius(Z=3, N=2)`` keeps neither ``Z`` nor
# ``N`` on the constructed object, so the loaded rate law cannot say where its
# numbers came from. That is written up as a database/engine design finding in
# ``docs/i134-sole-reaction-selection.md``.

#: Kelvin per electronvolt, from the CODATA-2018 elementary charge (exact by
#: definition) and RMG's Boltzmann constant. Physical constants are not what is
#: under test here; the tolerance below is wide enough to absorb the difference
#: between any two modern constant sets and still ~10^95 tighter than what is
#: needed to see an ``A`` replaced by ``1e100``.
KELVIN_PER_EV = 1.602176634e-19 / constants.kB

#: Relative tolerance for every anchor. Measured, not chosen: the worst residual
#: between these restatements and the shipped rate laws is 1.55e-7, on the Voronov
#: form, where the eV-to-K conversion differs between constant sets; the two other
#: fits agree to 5e-16, which is float noise. 1e-5 is a 64x margin over the worst
#: of those and ~10^95 short of the mutation this exists to catch.
#: ``evidence/ANCHOR-agreement.stdout.log`` is the measurement and
#: ``evidence/anchor_agreement.py`` re-runs it against whatever this file
#: currently asserts, so the margin cannot silently stop being true.
ANCHOR_RTOL = 1.0e-5


def _cm3_per_molecule_to_si(value):
    """cm^3 molecule^-1 s^-1 -> m^3 mol^-1 s^-1, the units RMG evaluates in."""
    return value * 1.0e-6 * constants.Na


def _badnell_rr(temperature, A, B, T0, T1, C, T2):
    """Badnell (2006), A&A 447, 389, eq. 1 -- the radiative-recombination fit.

        alpha(T) = A [ sqrt(T/T0) (1 + sqrt(T/T0))^(1-B')
                                  (1 + sqrt(T/T1))^(1+B') ]^-1
        B' = B + C exp(-T2/T)

    `A` is in cm^3 s^-1 and the temperatures in K, as the paper tabulates them.
    """
    b = B + C * math.exp(-T2 / temperature)
    root0 = math.sqrt(temperature / T0)
    root1 = math.sqrt(temperature / T1)
    return A / (root0 * (1.0 + root0) ** (1.0 - b) * (1.0 + root1) ** (1.0 + b))


def _voronov_ei(temperature, A, dE, P, K, X):
    """Voronov (1997), ADNDT 65, 1, eq. 1 -- the electron-impact ionisation fit.

        k(Te) = A (1 + P sqrt(U)) / (X + U) U^K exp(-U),   U = dE / kTe[eV]

    `A` is in cm^3 s^-1, `dE` in eV, `temperature` in K.
    """
    reduced = dE / (temperature / KELVIN_PER_EV)
    return (A * (1.0 + P * math.sqrt(reduced)) / (X + reduced)
            * reduced ** K * math.exp(-reduced))


def _shull_van_steenberg_rr(temperature, A_rad, eta):
    """Shull & Van Steenberg (1982), ApJS 48, 95, eq. 1 -- radiative recombination.

        alpha_r(T) = A_rad (T / 1e4 K)^-eta

    `A_rad` is in cm^3 s^-1. Row AR1 of their Table 2 is argon's.
    """
    return A_rad * (temperature / 1.0e4) ** -eta


class RateAnchor:
    """One published rate value the entry's own rate law has to reproduce.

    `expected` is in m^3/(mol*s) -- the SI units the reactor's evaluator returns
    -- and `source` is the sentence a reader of a red run needs in order to know
    which publication disagrees with the database.

    `temperature` is an ELECTRON temperature and `gas_temperature` is the GAS
    temperature the same evaluation is made at. Every published form restated in
    this file is a function of Te, and the entry's rate law is evaluated against
    it through :func:`evaluate_as_the_solver_does` at ``(gas_temperature, Te)``
    -- the reactor's own dispatch -- rather than through the one-temperature
    interface, which on ``TwoTemperaturePlasma`` would hide every gas-temperature
    term the entry carries and the published form does not.

    **Why the gas temperature is a field and not a constant.** Until this round
    every anchor was taken at :data:`SOLVER_TGAS`, so the whole anchor set lay on
    a line of constant ``Tgas`` -- and on that line ``TwoTemperaturePlasma``'s
    ``exp(-Ea_g/(R*T))`` is a constant that is absorbed into ``A`` exactly. The
    two parameters were not independently determined by any number of anchors,
    only by the number of *distinct gas temperatures*, which was one. Measured on
    the built module, against a mutation that sets ``Ea_g = 10 kJ/mol`` and
    divides ``A`` by ``exp(-Ea_g/(R*298.15))``:

        Tgas=298.15  Te=2e4   good=1.445845e+05  mut=1.445845e+05  ratio 1.000000000
        Tgas=298.15  Te=2e5   good=3.229400e+04  mut=3.229400e+04  ratio 1.000000000
        Tgas=1000    Te=2e4   good=1.445845e+05  mut=2.453061e+06  ratio 16.966
        Tgas=5000    Te=2e4   good=1.445845e+05  mut=6.420619e+06  ratio 44.407

    Identical to nine decimals where the anchors looked, and 17x to 44x wrong one
    axis over. :data:`ANCHOR_GAS_TEMPERATURES` and
    :func:`anchor_sensitivity` are the repair and its proof respectively.
    """

    def __init__(self, temperature, expected, source,
                 gas_temperature=SOLVER_TGAS, rtol=ANCHOR_RTOL):
        self.temperature = temperature
        self.gas_temperature = gas_temperature
        self.expected = expected
        self.source = source
        self.rtol = rtol

    def __str__(self):
        return '(Tgas = {0} K, Te = {1} K)'.format(self.gas_temperature,
                                                   self.temperature)

    def agrees(self, found):
        """True when `found` matches, and False when it is NaN.

        The comparison is ``<=``, whose failure branch is the one NaN takes, so a
        rate law that evaluates to NaN is REPORTED rather than passed. That is the
        opposite polarity from the validity-window comparison in
        :meth:`ShippedEntry.assert_as_shipped`, which had to be repaired this
        round -- see :func:`nan_polarity_census`.
        """
        return abs(found - self.expected) <= self.rtol * abs(self.expected)


#: Badnell (2006) row for the Li II -> Li I stage, (Z=3, N=2), as the entry's own
#: longDesc quotes it from the paper.
_BADNELL_LITHIUM = dict(A=8.7e-12, B=0.364, T0=147.0, T1=7.153e6,
                        C=0.1508, T2=7.154e5)
#: Voronov (1997) row for neutral lithium, (Z=3, N=3).
_VORONOV_LITHIUM = dict(A=1.39e-7, dE=5.4, P=0.0, K=0.41, X=0.438)
#: Shull & Van Steenberg (1982) Table 2 row AR1, argon's radiative recombination.
_SVS_ARGON = dict(A_rad=3.77e-13, eta=0.651)

#: The Boltzmann constant in eV/K that ``VoronovEIArrhenius.populate_from_yaml``
#: uses to convert that table's eV validity bounds into K. It is NOT the constant
#: :data:`KELVIN_PER_EV` is built from -- the two differ in the 10th digit -- so
#: the Voronov window below is stated through this one, which is the one the
#: shipped object's ``Tmin``/``Tmax`` actually came out of.
_VORONOV_KB_EV_PER_K = 8.617333262e-5

#: SIX electron temperatures, and the count is derived rather than chosen -- see
#: :func:`anchor_grid` for what makes an anchor set adequate.
#:
#: The binding constraint is the widest rate law. ``BadnellRRArrhenius`` has six
#: free parameters, and an anchor set can separate at most as many independent
#: directions as it has points. Two (what this file had) determined nothing beyond
#: a magnitude and a slope; measured, three give the Badnell fit rank 3 of 6, and
#: six give it rank 6 of 6.
#:
#: All six lie inside every entry's own declared validity window, whose
#: intersection is [11604.5 K, 1e7 K]: Voronov's Tmin is 11604.5 K (which is why
#: 2e4 K and not 1e4 K is the low anchor) and Badnell's Tmax is 1e7 K (which is
#: why the high anchor is 5e6 K and not higher). They are spread roughly
#: geometrically across that band, because the parameters they have to separate
#: enter as powers and exponentials of Te rather than linearly.
ANCHOR_TEMPERATURES = (2.0e4, 6.0e4, 2.0e5, 6.0e5, 2.0e6, 5.0e6)

#: TWO gas temperatures, for the entries whose rate law is a function of one.
#:
#: 298.15 K is :data:`SOLVER_TGAS`, the 5 torr argon deck's own gas temperature.
#: 2000 K is a second one a low-pressure argon discharge plausibly reaches, and
#: its exact value is not load-bearing: measured, the smallest non-degenerate
#: singular value of the argon design matrix is flat to five digits for a second
#: gas temperature anywhere in 1000-5000 K. What IS load-bearing is that there
#: are two of them, because with one the ``exp(-Ea_g/(R*T))`` factor is a
#: constant absorbed into ``A`` exactly.
#:
#: This axis is added to an entry's anchors only when the entry's rate law
#: actually varies with the gas temperature, and that is MEASURED on the rate law
#: rather than keyed on its class name -- see :func:`rate_law_admits_gas_temperature_dependence`.
#: For the two lithium entries, whose published forms and whose implementations
#: are both functions of Te alone, a second gas temperature would restate the
#: first anchor to the last bit and assert nothing; measured, Badnell and Voronov
#: both return bit-identical k at (298.15, 2e4) and at (3000, 2e4).
ANCHOR_GAS_TEMPERATURES = (SOLVER_TGAS, 2000.0)

#: The relative step used to differentiate a rate law with respect to one of its
#: own parameters in :func:`anchor_sensitivity`. Central differences, so the
#: truncation error is O(h^2) and the rounding error O(eps/h); 1e-6 sits near the
#: minimum of that sum for a smooth double-precision function, and the measured
#: sensitivities below are stable to five digits across h = 1e-5 .. 1e-7.
SENSITIVITY_STEP = 1.0e-6


def _varies_with_gas_temperature(kinetics, rtol=1.0e-12):
    """Whether THIS rate law object returns different k at two gas temperatures."""
    low, high = ANCHOR_GAS_TEMPERATURES
    probe = ANCHOR_TEMPERATURES[0]
    at_low = evaluate_as_the_solver_does(kinetics, probe, low)
    at_high = evaluate_as_the_solver_does(kinetics, probe, high)
    return abs(at_high - at_low) > rtol * abs(at_low)


def rate_law_admits_gas_temperature_dependence(kinetics, step=1.0e-3):
    """Whether `kinetics`'s FORM has a gas-temperature parameter -- not whether
    this particular instance happens to use it.

    **The distinction is the whole point, and getting it wrong was a real defect
    in the first draft of this round's repair.** The obvious implementation asks
    the object in hand whether ``k`` changes between two gas temperatures. On the
    shipped argon entry that returns ``False``, because the entry declares
    ``Ea_g = Ea_e = 0`` and is therefore genuinely gas-temperature invariant. A
    grid built on that answer would drop the second gas-temperature axis from the
    one entry that needs it -- and it would drop it precisely BECAUSE the entry is
    currently correct, so the axis would be absent exactly when a mutation that
    sets ``Ea_g = 10 kJ/mol`` arrived. The discriminator would remove the check
    that catches the bad state on the evidence that the good state does not need
    it.

    That is a defect class worth naming, because it is not specific to
    temperatures: **a discriminator evaluated on the good state reports the
    property the good state has, and switches off the check that would have seen
    the bad one.** The test for it is to ask whether the bad state is REACHABLE,
    not whether the current state is bad.

    So the question asked here is whether ANY perturbation of the rate law's own
    parameters can introduce a gas-temperature dependence. If one can, the
    parameter exists in the functional form, the entry's invariance is a claim
    rather than a tautology, and the anchors need a second gas temperature both
    to check that claim and to separate ``A`` from ``Ea_g``. Measured: true for
    ``TwoTemperaturePlasma`` (perturbing ``Ea_g`` or ``Ea_e`` does it), false for
    ``BadnellRRArrhenius`` and ``VoronovEIArrhenius``, whose forms have no gas
    temperature in them at all and for which no perturbation of any parameter can
    produce one.

    Still not keyed on the class name, for the reason
    :func:`evaluate_as_the_solver_does` gives: a check that said
    ``isinstance(kinetics, TwoTemperaturePlasma)`` would go quietly blind on the
    next two-temperature class the database learns to author.
    """
    if _varies_with_gas_temperature(kinetics):
        return True
    for name in rate_law_parameter_names(kinetics):
        for direction in (step, -step):
            try:
                moved = perturb_rate_law_parameter(kinetics, name, direction)
            except Exception:               # pragma: no cover - a refused parameter
                continue
            if _varies_with_gas_temperature(moved):
                return True
    return False


def anchor_grid(kinetics):
    """The ``(Tgas, Te)`` points an entry carrying `kinetics` is anchored at.

    Three electron temperatures always; two gas temperatures when -- and only
    when -- the rate law is measurably a function of the gas temperature. Six
    points for the argon entry, three for each lithium entry.

    **Why this shape, and why the obvious fix is not it.** Take logs of the
    Kossyi form production evaluates for ``TwoTemperaturePlasma``::

        k    = A * (Te/T0)^n * exp(-Ea_g/(R*T)) * exp(Ea_e*(Te-T)/(R*T*Te))
        ln k = ln A - n*ln T0 + n*ln Te + (Ea_e - Ea_g)/(R*T) - Ea_e/(R*Te)

    ``ln k`` is LINEAR in the unknowns over the basis ``{1, ln Te, 1/T, 1/Te}``,
    so an anchor set determines the parameters exactly when the design matrix
    whose rows are ``[1, ln Te, 1/T, -1/Te]`` has full column rank. Measured on
    the built module:

        ===========================================  ======  ====  ==========
        design (rank over A, n, Ea_g, Ea_e)          points  rank  sigma_min
        ===========================================  ======  ====  ==========
        1 Tgas x 2 Te   (what this file had)              2     2  --
        2 Tgas x 2 Te   (the obvious fix)                 4     3  8.7e-18
        3 Te at one Tgas + 1 point at a second            4     4  1.4e-3
        3 Te x 2 Tgas                                     6     4  2.0e-3
        6 Te x 2 Tgas   (this design)                    12     4  2.9e-3
        ===========================================  ======  ====  ==========

    Six electron temperatures rather than three because the argon entry is not
    the widest rate law here -- Badnell's has six parameters and needs six points
    on the Te axis alone. See :data:`ANCHOR_TEMPERATURES`.

    **The 2x2 grid is the trap, and it is worth naming.** The model is additive
    in ``ln Te`` and ``1/T`` with no interaction term between them, so a tensor
    product of two values on each axis is rank-deficient however far apart the
    two values are: its fourth row is identically the second plus the third minus
    the first. Adequacy needs at least three distinct ``Te`` AND at least two
    distinct ``Tgas``. "Add another point" is not a design, and adding points one
    at a time is the process that produced both of the last two defects here.

    Rank is necessary and not sufficient: a formally full-rank design with a bad
    condition number still hides a perturbation below :data:`ANCHOR_RTOL`. What
    this suite therefore asserts is not a count of points but the measured
    sensitivity of the anchor set to each parameter of each entry's own rate law
    -- see :func:`anchor_sensitivity` and
    ``TestTheAnchorSetDeterminesEveryParameter``.
    """
    if rate_law_admits_gas_temperature_dependence(kinetics):
        return tuple((gas, electron)
                     for gas in ANCHOR_GAS_TEMPERATURES
                     for electron in ANCHOR_TEMPERATURES)
    return tuple((SOLVER_TGAS, electron) for electron in ANCHOR_TEMPERATURES)


class _UnusedSlot:
    """A ``__reduce__`` slot the class hard-codes as a literal ``None``.

    A distinct object rather than ``None`` itself, so that "this file records this
    slot as unused" and "the census forgot to say what this slot is" cannot be
    confused, and so that no attribute name can ever collide with it.
    """

    def __repr__(self):
        return 'UNUSED_SLOT'


UNUSED_SLOT = _UnusedSlot()

#: Sentinel for ``getattr`` -- distinct from ``None``, which several of these
#: fields legitimately hold.
_NO_SUCH_ATTRIBUTE = object()

#: The free parameters of each shipped plasma rate law: the name each one is
#: readable under, its position in the class's own ``__reduce__`` argument tuple,
#: and what every other slot of that tuple carries.
#:
#: **Why the position and the arity are both here.** A rate law's parameters have
#: to be perturbed generically -- one at a time, without a hand-written mutator
#: per class -- and ``inspect.signature`` cannot see a compiled ``cdef class``
#: (measured: it reports ``(self, *args, **kwargs)`` for all three). ``__reduce__``
#: is the one interface that does enumerate them, because pickling has to
#: round-trip every field that matters. Rebuilding ``cls(*args)`` with one slot
#: replaced is therefore a perturbation of exactly one parameter, expressed in the
#: class's own terms.
#:
#: This table is a hand-written census of another module's fields -- the mirror
#: shape this campaign keeps finding stale -- so it is tied back by
#: ``test_the_parameter_census_matches_each_rate_laws_own_arity``.
#:
#: **Every slot is accounted for, and the recorded arity is gone.** Until this
#: round the tie-back was the total arity plus "each recorded parameter slot is in
#: range and readable". Both survive the case that matters: a **bookkeeping slot
#: replaced by a new functional parameter** keeps the arity the same and leaves
#: every recorded parameter slot valid, so a new free parameter could appear in
#: one of these rate laws and escape every assertion in this file -- which is the
#: precise failure the census was written to prevent, one level down.
#:
#: So ``reduce_slots`` names **all** of them, in ``__reduce__`` order: an
#: attribute name for a slot that carries a field, and :data:`UNUSED_SLOT` for one
#: the class hard-codes as a literal ``None`` (Badnell's three reserved slots and
#: Voronov's ``(Z, N, yaml_path_or_obj)`` lookup triple). The arity is now
#: ``len(reduce_slots)`` -- derived, not a second hand-written number that could
#: disagree with the first. :func:`assert_reduce_slots_match` checks each slot
#: against the value read off the live object under the recorded name, so a slot
#: whose meaning changed is red even where the count did not.
#:
#: The non-parameter slots are the bookkeeping fields (``electrons``,
#: ``Tmin``/``Tmax``/``Pmin``/``Pmax``, ``uncertainty``, ``solute``, ``comment``),
#: asserted by name elsewhere in this file and not part of the functional form.
RATE_LAW_PARAMETERS = {
    'TwoTemperaturePlasma': dict(
        reduce_slots=('A', 'n', 'Ea_g', 'Ea_e', 'T0',
                      'electrons',
                      'Tmin', 'Tmax', 'Pmin', 'Pmax',
                      'uncertainty', 'solute', 'comment'),
        parameters=(('A', 0), ('n', 1), ('Ea_g', 2), ('Ea_e', 3), ('T0', 4))),
    'BadnellRRArrhenius': dict(
        reduce_slots=('A', 'B', 'T0', 'T1', 'C', 'T2',
                      'electrons',
                      UNUSED_SLOT, UNUSED_SLOT, UNUSED_SLOT,
                      'Tmin', 'Tmax', 'Pmin', 'Pmax',
                      'uncertainty', 'solute', 'comment'),
        parameters=(('A', 0), ('B', 1), ('T0', 2), ('T1', 3), ('C', 4), ('T2', 5))),
    'VoronovEIArrhenius': dict(
        reduce_slots=('A', 'P', 'X', 'K', 'dE',
                      'electrons',
                      UNUSED_SLOT, UNUSED_SLOT, UNUSED_SLOT,
                      'Tmin', 'Tmax', 'Pmin', 'Pmax',
                      'uncertainty', 'solute', 'comment'),
        parameters=(('A', 0), ('P', 1), ('X', 2), ('K', 3), ('dE', 4))),
}


def assert_reduce_slots_match(kinetics, census=None):
    """Every ``__reduce__`` slot of `kinetics` is the field this file records there.

    The tie-back for :data:`RATE_LAW_PARAMETERS`. Three claims, and the second is
    the one that is new:

    1. the tuple has as many slots as the census names -- the old arity tripwire,
       now derived from the slot map rather than carried as a second number;
    2. **every slot holds what the census says it holds**, checked by reading the
       named attribute off the live object and comparing it to the value at that
       position. A slot recorded as :data:`UNUSED_SLOT` must still be a literal
       ``None``. This is what makes a *changed* slot visible: a bookkeeping slot
       turned into a functional parameter leaves the arity alone and leaves every
       recorded parameter slot valid, and it stops being ``None``;
    3. the ``parameters`` list agrees with the slot map about where each free
       parameter lives, so the two halves of the census cannot drift apart --
       ``perturb_rate_law_parameter`` rebuilds ``cls(*args)`` from that list, and a
       parameter pointed at the wrong slot would perturb the wrong field and every
       sensitivity measured from it would be a measurement of something else.

    Comparison is ``is`` first and ``repr`` as a fallback. Identity holds for every
    slot of all three shipped rate laws except Voronov's ``dE``, which is carried
    in the tuple as ``self._dE_eV`` -- a bare float, so the property returns an
    equal object rather than the same one. Measured, not assumed.

    :param census: the census row to check against; defaults to the one recorded
        for this class. Passing one explicitly is how
        ``test_a_bookkeeping_slot_that_became_a_parameter_is_refused`` puts a
        stand-in class through this same function instead of a copy of it.
    """
    name = type(kinetics).__name__
    census = RATE_LAW_PARAMETERS[name] if census is None else census
    slots = census['reduce_slots']
    arguments = kinetics.__reduce__()[1]
    assert len(arguments) == len(slots), (
        '{0} now reduces to {1} arguments and this file accounts for {2}. Something '
        'was added to or removed from the rate law; work out whether it is a free '
        'parameter of the functional form (then add it to RATE_LAW_PARAMETERS with '
        'its slot, and give it a pin) or bookkeeping (then only the slot map '
        'moves).'.format(name, len(arguments), len(slots)))

    for position, descriptor in enumerate(slots):
        held = arguments[position]
        if descriptor is UNUSED_SLOT:
            assert held is None, (
                '{0}.__reduce__() slot {1} carries {2!r}, and this file records that '
                'slot as one the class hard-codes as None. A bookkeeping slot that '
                'has become a field keeps the ARITY unchanged and leaves every '
                'recorded parameter slot in range and readable, so neither of those '
                'checks can see it -- which is why every slot is accounted for here '
                'and not only the parameters. Work out what now lives there: if it '
                'is a free parameter of the functional form it needs a name in '
                'reduce_slots, an entry in parameters, and a pin.'.format(
                    name, position, held))
            continue
        attribute = getattr(kinetics, descriptor, _NO_SUCH_ATTRIBUTE)
        assert attribute is not _NO_SUCH_ATTRIBUTE, (
            '{0} no longer exposes {1!r}, so the census names a field at slot {2} '
            'that cannot be read off the object at all.'.format(name, descriptor, position))
        assert held is attribute or repr(held) == repr(attribute), (
            '{0}.__reduce__() slot {1} is recorded as {2!r}, but that slot holds '
            '{3!r} while {2} reads {4!r}. The slot map has drifted from the class: '
            'either the fields were reordered, or this slot now carries something '
            'else.'.format(name, position, descriptor, held, attribute))

    for parameter, position in census['parameters']:
        assert 0 <= position < len(slots), (
            '{0}\'s recorded slot {1} for {2} is outside its reduce tuple'.format(
                name, position, parameter))
        assert slots[position] == parameter, (
            '{0}\'s parameter list puts {1!r} at slot {2}, and the slot map says that '
            'slot holds {3!r}. perturb_rate_law_parameter() rebuilds the class from '
            'that slot, so the two disagreeing means every sensitivity measured for '
            '{1} was measured on a different field.'.format(
                name, parameter, position, slots[position]))


def rate_law_parameter_names(kinetics):
    """The names of `kinetics`'s own rate-law parameters, in ``__reduce__`` order."""
    census = RATE_LAW_PARAMETERS[type(kinetics).__name__]
    return tuple(name for name, _ in census['parameters'])


def perturb_rate_law_parameter(kinetics, name, delta):
    """A copy of `kinetics` with ONE parameter moved by `delta`, and nothing else.

    `delta` is dimensionless and is applied against a per-parameter SCALE: the
    magnitude of the parameter's own authored value where that is non-zero, and
    1.0 in the parameter's own authored units where it is zero. So ``delta`` reads
    as a relative perturbation for ``A``, ``n``, ``T0``, ``T1``, ``T2``, ``B``,
    ``C``, ``K``, ``X`` and ``dE``, and as an absolute perturbation in kJ/mol for
    the argon entry's ``Ea_g`` and ``Ea_e`` (both authored as 0.0, where a
    relative step would be a no-op) and in bare units for Voronov's ``P`` (0.0).
    That rule is stated rather than chosen per parameter, so a parameter that
    later becomes non-zero does not silently change what a reported threshold
    means.

    Rebuilt through ``__reduce__`` rather than by assigning to the attribute:
    assignment goes through each class's property setter, which re-derives
    dependent state in ways that differ between the three classes, and a
    perturbation that also moved something else would make every sensitivity
    below a measurement of the wrong thing.
    """
    census = RATE_LAW_PARAMETERS[type(kinetics).__name__]
    slots = dict(census['parameters'])
    if name not in slots:
        raise KeyError('{0} has no rate-law parameter {1!r}; it has {2}'.format(
            type(kinetics).__name__, name, ', '.join(slots)))
    cls, arguments = kinetics.__reduce__()
    arguments = list(arguments)
    held = arguments[slots[name]]
    if isinstance(held, float):             # Voronov's dE, carried as a bare float
        value, units = held, None
    else:
        value, units = held.value, held.units
    scale = abs(value) if value != 0.0 else 1.0
    moved = value + delta * scale
    arguments[slots[name]] = moved if units is None else (moved, units)
    return cls(*arguments)


def anchor_sensitivity(kinetics, grid=None, step=SENSITIVITY_STEP):
    """How sharply an anchor set sees each parameter of `kinetics`'s rate law.

    Returns ``(names, points, jacobian)``, where ``jacobian[i][j]`` is
    ``d ln k / d delta_j`` at anchor point ``i`` -- the fractional change in the
    rate coefficient per unit perturbation of parameter ``j`` in the units
    :func:`perturb_rate_law_parameter` defines. Central differences.

    **Why log-sensitivity and not the rate itself.** :meth:`RateAnchor.agrees`
    compares ``abs(found - expected) <= rtol * abs(expected)``, i.e. it tests
    ``abs(k_found/k_expected - 1) <= rtol``, and to first order that quantity IS
    ``abs(d ln k)``. So a column of this matrix converts directly into the
    smallest perturbation of that parameter the anchor set refuses, with no
    second constant to choose:

        smallest detected delta_j  =  ANCHOR_RTOL / max_i abs(jacobian[i][j])

    and a parameter whose column is identically zero is one NO anchor set of this
    shape can ever see, at any tolerance and any number of points.

    **Why a numerical Jacobian rather than the closed form.** Only
    ``TwoTemperaturePlasma`` is linear in its parameters after taking logs.
    ``BadnellRRArrhenius`` and ``VoronovEIArrhenius`` are not, so there is no
    exact design matrix to write down for them, and hand-differentiating two
    published fits inside a test file would be a third transcription to keep in
    step. Differentiating the built objects through production's own evaluator
    keeps the count of things that can drift at one. The cost is a step size, and
    :data:`SENSITIVITY_STEP` records how it was chosen and over what range the
    answer is stable.
    """
    if grid is None:
        grid = anchor_grid(kinetics)
    names = rate_law_parameter_names(kinetics)
    jacobian = np.zeros((len(grid), len(names)))
    for column, name in enumerate(names):
        up = perturb_rate_law_parameter(kinetics, name, step)
        down = perturb_rate_law_parameter(kinetics, name, -step)
        for row, (gas, electron) in enumerate(grid):
            jacobian[row, column] = (
                math.log(evaluate_as_the_solver_does(up, electron, gas))
                - math.log(evaluate_as_the_solver_does(down, electron, gas))
            ) / (2.0 * step)
    return names, tuple(grid), jacobian


def smallest_detected_perturbations(kinetics, grid=None, rtol=ANCHOR_RTOL):
    """Per parameter, the smallest lone perturbation the anchors refuse.

    ``None`` where the anchors can never see the parameter at all, however large
    the perturbation -- which for the argon entry's ``T0`` is the true answer and
    not a tolerance to be tightened. See
    ``test_a_parameter_no_anchor_can_reach_is_pinned_by_a_direct_assertion``.
    """
    names, _, jacobian = anchor_sensitivity(kinetics, grid=grid)
    detected = {}
    for column, name in enumerate(names):
        reach = float(np.max(np.abs(jacobian[:, column])))
        detected[name] = (rtol / reach) if reach > 0.0 else None
    return detected


def _matches_within(found, expected, rtol):
    """``found`` equals ``expected`` to `rtol`, with NaN and inf reported as wrong.

    **The defect this exists to stop, measured.** The validity-window check used
    to read::

        elif found is None or abs(found - expected) > 1e-9 * abs(expected):

    -- a comparison whose FAILURE branch is the one NaN takes. ``abs(nan - 1e4) >
    1e-9 * 1e4`` is ``False``, so a rate law declaring ``Tmin = nan`` passed the
    window assertion, and it was not caught anywhere else either: measured,
    ``max(1e4, nan)`` returns ``1e4`` and ``min(1e8, nan)`` returns ``1e8``, so
    :func:`electron_temperature_response` clamps to the unchanged band and keeps
    its answer. Meanwhile production's own ``is_temperature_valid(20000)`` returns
    ``False`` on that same object -- the entry is unusable and every check here
    said it was fine.

    So the polarity is inverted deliberately: this asks whether the value IS
    right, and anything that is not a finite number in range is not. A comparison
    phrased as "is it wrong" hands NaN a free pass; one phrased as "is it right"
    does not.
    """
    if found is None:
        return False
    if not math.isfinite(found):
        return False
    if expected is None:
        return False
    if not math.isfinite(expected):             # pragma: no cover - no such entry
        return False
    return abs(found - expected) <= rtol * abs(expected)


#: Every numeric comparison in this file that decides whether a field is accepted,
#: and what each one does when handed NaN. Recorded because the defect repaired
#: this round was ONE comparison of the wrong polarity, and the class is broader
#: than the line: a check phrased as "is this value wrong" passes NaN, and one
#: phrased as "is this value right" refuses it.
#:
#: This table is a claim about behaviour, so it is tied back rather than left as
#: prose: ``test_the_nan_polarity_census_is_true_of_the_real_comparisons`` drives
#: each of these comparisons with NaN and asserts the recorded verdict. A
#: comparison that changes polarity turns that check red.
NAN_POLARITY_CENSUS = {
    'RateAnchor.agrees': (
        'REFUSES',
        'abs(found - expected) <= rtol * abs(expected) -- an "is it right" '
        'comparison, so NaN takes the False branch and the anchor is reported.'),
    'ShippedEntry.rate_law_parameters via _matches_within': (
        'REFUSES',
        'delegates to _matches_within, which tests math.isfinite explicitly '
        'before comparing.'),
    'ShippedEntry.temperature_window via _matches_within': (
        'REFUSES',
        'was `abs(found - expected) > 1e-9 * abs(expected)`, an "is it wrong" '
        'comparison that NaN passed; now _matches_within, which refuses it. '
        'This is the defect this round repaired.'),
    'electron_temperature_response band clamp': (
        'PASSES',
        'max(TE_BAND[0], Tmin) and min(TE_BAND[1], Tmax) both return the '
        'TE_BAND end when the other operand is NaN, so the response is read '
        'over the unchanged band and comes back unchanged. This is why the '
        'window assertion -- and not the response -- has to be the check that '
        'sees a NaN bound.'),
}


class ShippedEntry:
    """Which entry a check is about, and what that entry has to carry.

    Two halves, and keeping them apart is the whole point of this class.

    **The provenance half -- what SELECTS.** ``library`` is a library label, bound
    to a source file by :class:`LoadedLibraries`, and ``entry_index`` is the
    ``index =`` field authored on one entry in it. Nothing else. That pair is the
    entry's identity in the database, and it is complete by construction: an entry
    cannot be in two libraries, and within one library
    ``KineticsLibrary.load_entry`` refuses to load a second entry under an index
    already taken (``assert index not in self.entries``,
    ``rmgpy/data/kinetics/library.py:680``), so the index is unique *by the
    loader's own enforcement* and ``library.entries`` is keyed on it.

    **Why the index and not the label -- the correction this round makes.** The
    previous handle was ``entry.label``, and in this database format an entry's
    label **is its reaction string**: ``KineticsLibrary.load`` parses
    ``entry.label`` for the participants, for the ``(+M)`` collider and for the
    direction, and refuses any entry whose ``reversible`` field disagrees with the
    arrow in its own label (``library.py:571``). So the label is not a name for the
    entry at all -- it is a *rendering of the entry's chemistry*, and it moves
    whenever the chemistry moves. That made the handle inherit every mutability the
    assertions exist to catch: a wrong-chemistry entry failed as a **selection
    miss** ("no entry labelled ...") instead of as an assertion naming what was
    wrong, and a later entry taking the vacated label could redirect the selection
    outright.

    Measured on this database, in ``evidence/HANDLE-invariance.stdout.log``: there
    is **no chemistry-preserving rename** of these entries, because every loadable
    label is a valid reaction string for the chemistry the entry declares. Rename
    and chemistry-change are the same operation in this format. The index, by
    contrast, is an authored integer that renders nothing about the chemistry, and
    the same measurement shows it unmoved by a rename, a reorder of the entries in
    ``reactions.py`` and a reversibility flip alike.

    **The physics half -- what is ASSERTED.** ``entry_label``, ``order``,
    ``response``, ``reversible``, ``electrons``, ``kinetics_class``,
    ``shipped_fit``, ``rate_anchors`` and the three authorable flags are what the
    selected entry must then turn out to carry. :meth:`assert_as_shipped` checks
    each one separately and names the attribute that is wrong. They are not part of
    the key, and putting any of them back into it would reintroduce exactly the
    defect this class exists to remove: a wrong value would once again pick a
    different entry instead of failing.

    **The label is now on this side of the line**, which is the whole point of the
    move: an entry whose chemistry has been rewritten is still *that entry*, is
    still selected, and is refused by an assertion that names the label, the
    reversibility and the rate together.

    **Why the inversion was needed.** The previous key was the participants'
    graphs, the rate order and the sense of the Te response, and every round
    found a substitution it accepted -- because the participant lists do not
    carry the electron, so two different channels over the same heavy species are
    written identically in this database::

        radiative recombination     Li+ + e-   =>  Li + hv    second order
        three-body recombination    Li+ + 2 e- =>  Li + e-    third order

    Both are stored as ``[Lip] => [Li]`` with ``electrons = -1``, both are owned
    by ``PlasmaRadiativeRecombination`` and so carry the same ``(1, 0)``
    placement declaration. Adding the order separated those two; then the shipped
    *ionisation* fit ``VoronovEIArrhenius(Z=3, N=3, electrons=-1)`` pasted onto
    the recombination entry got through, because it is second order too. Adding
    the Te response separated that one; then a ``reversible = True`` radiative
    entry, a generic order-2 falling ``ElectronCollisionPlasma``, and the shipped
    *argon* ``TwoTemperaturePlasma`` fit pasted onto lithium all got through. The
    attribute space is open. For any finite attribute key, a rate law that agrees
    on every attribute in it can be written, so the attribute route cannot
    converge -- and each round of widening made one more physical fact into a
    silent selector instead of a loud assertion.

    All five of those substitutions are now red, and each is red as an assertion
    naming the attribute it got wrong. The evidence is in
    ``docs/i134-sole-reaction-selection/evidence/BYPASS-*``.

    **What a database-side rename, reorder or renumber does to this key.**
    Measured, not reasoned about; ``evidence/HANDLE-invariance.stdout.log`` is the
    log and ``evidence/handle_invariance.py`` re-runs it.

    * A **rename** -- which in this format is always also a chemistry change --
      leaves the selection exactly where it was, and the entry is refused by
      assertion, naming the label and whatever else moved with it. Under the old
      label handle this was a selection miss that never reached an assertion.
    * A **reorder** of the entries in ``reactions.py`` does nothing: position is
      not part of the key, which is the original ``reactions[0]`` defect closed
      for good.
    * A **renumber** moves the handle, and is therefore the one edit that makes
      this key fail. It fails *loudly*: the index is absent and
      :func:`select_entry` raises, naming every index and label the library does
      hold. That is the same shape of failure a rename used to produce under the
      old handle, now confined to an edit that changes no chemistry -- which is
      the trade this round makes deliberately, because the edit that changes no
      chemistry is the rare one and the edit that changes chemistry is the one
      the assertions exist for.
    * A **swap of two indices** is the residual risk, and it is the only edit that
      can silently *redirect* the selection. It does not stay silent: the
      redirected entry is a different reaction, so the reactants, the products and
      the fit all go red by name -- measured, four named attributes for the
      lithium/argon swap in this library. The conversion of that redirect into a
      loud assertion is exactly why ``entry_label`` moved to the asserted half.

    **The collision the label handle admitted is gone.** Two entries of one library
    may legitimately share a label -- a duplicate-marked pair -- and the old handle
    had to refuse in that case rather than pick. Two entries cannot share an index:
    the loader refuses the file. ``test_no_shipped_plasma_entry_label_collides``
    stays, because a label collision is still a fact worth knowing about the
    database, but it no longer bounds what this key can name.
    """

    def __init__(self, name, library, entry_index, entry_label, reactants, products,
                 order, response, reversible, electrons,
                 kinetics_class, temperature_window, shipped_fit=None, fit_source=None,
                 rate_anchors=(), rate_law_parameters=None, anchor_degeneracies=(),
                 duplicate=False, allow_pdep_route=False, elementary_high_p=False,
                 allow_max_rate_violation=False):
        self.name = name
        #: PROVENANCE. The label of the library the entry lives in.
        self.library = library
        #: PROVENANCE. The authored ``index =`` field of the entry within that
        #: library, which is also its ``library.entries`` key.
        self.entry_index = entry_index
        #: Everything below is ASSERTED on the selected entry, never used to
        #: select it -- including the entry's own label, which in this database
        #: format renders the chemistry and so belongs here and not in the key.
        self.entry_label = entry_label
        self.reactants = tuple(reactants)
        self.products = tuple(products)
        self.order = order
        self.response = response
        self.reversible = reversible
        self.electrons = electrons
        self.kinetics_class = kinetics_class
        #: ASSERTED. ``(Tmin, Tmax)`` in K as the shipped rate law declares them.
        #: Required, not defaulted, and not ``None``-able: an entry whose window is
        #: genuinely open would have to say so as ``(None, None)`` deliberately.
        #:
        #: **Why it needs its own assertion, when the Te response already reads
        #: these two fields.** It reads them to CLAMP to them
        #: (:func:`electron_temperature_response` intersects :data:`TE_BAND` with
        #: the declared range), so narrowing the window moves the two Te values the
        #: response is read at and the response is unchanged -- the check adapts to
        #: the very edit it would have to catch. Measured on the argon entry:
        #: narrowing ``Tmin`` from 1e4 K to 9e5 K leaves the response ``FALLS`` and
        #: every other check green. ``is_identical_to`` does compare ``Tmin``/
        #: ``Tmax``, so ``shipped_fit`` sees a narrowing on the two lithium entries
        #: -- and argon has no ``shipped_fit``, which is where the hole was.
        #:
        #: A validity range is a scientific claim about where a fit may be used,
        #: and on the argon entry it is a claim its own ``longDesc`` is careful to
        #: qualify (grid membership, not an accuracy bound). Changing it is a
        #: decision, so it is red once, at the moment it is made.
        self.temperature_window = tuple(temperature_window)
        #: The three flags an entry may author that nothing else here reads.
        #: All three are ``False`` on every entry this file names, measured; they
        #: are asserted so that setting them is a red check rather than a silent
        #: change in how RMG treats the reaction downstream.
        self.duplicate = duplicate
        self.allow_pdep_route = allow_pdep_route
        self.elementary_high_p = elementary_high_p
        #: The fourth authorable flag. Unlike the three above, this one is
        #: asserted on ``entry.item`` (the raw database entry), NOT on the
        #: constructed ``reaction`` the other three read: ``get_library_reactions``
        #: (``rmgpy/data/kinetics/library.py``, the ``else:  # pdep or standard
        #: library reaction`` branch) builds its ``LibraryReaction(...)`` without
        #: forwarding ``entry.item.allow_max_rate_violation``, so the constructed
        #: reaction always reads the ``False`` default regardless of what was
        #: authored. Asserting it on ``reaction`` would therefore be a check that
        #: can never fail -- exactly the defect class this file exists to catch.
        #: Measured ``False`` on every entry this file names.
        self.allow_max_rate_violation = allow_max_rate_violation
        #: MAGNITUDE. A tuple of :class:`RateAnchor`, each evaluating the entry's
        #: own rate law at one temperature and comparing it against a value
        #: computed here from the **published** closed form, with the published
        #: parameters transcribed into this file. See :class:`RateAnchor` for why
        #: this is the check ``shipped_fit`` cannot be.
        self.rate_anchors = tuple(rate_anchors)
        #: A zero-argument callable returning the rate law this entry is supposed
        #: to carry, constructed independently of the entry -- for the two lithium
        #: channels, straight out of the database's own ``badnell.yaml`` /
        #: ``voronov.yaml`` tables by ``(Z, N)``. That is the assertion that
        #: catches a fit for the *wrong element*, which no dimensional or
        #: monotonicity check can see. ``None`` where the shipped rate is authored
        #: in ``reactions.py`` rather than derived from a shipped table, because
        #: then the only available reference would be a hand-copied second copy of
        #: the database's own numbers, which drifts.
        self.shipped_fit = shipped_fit
        self.fit_source = fit_source
        #: ASSERTED. Each free parameter of the entry's own rate law, mapped to
        #: the value the PUBLICATION gives for it, in SI. Keys are the names in
        #: :data:`RATE_LAW_PARAMETERS`; the values are read off the object with
        #: ``value_si`` and compared exactly (to float noise).
        #:
        #: **Why this exists when the anchors already evaluate the rate law.**
        #: The anchors and this assertion see different things, and neither
        #: subsumes the other:
        #:
        #: * The anchors check the FORM and the DISPATCH -- that these fields,
        #:   put through production's own two-temperature evaluator, reproduce
        #:   the published number. A field assertion is blind to a rate law that
        #:   carries the right parameters and computes the wrong function of
        #:   them, or that production never routes to the right evaluator.
        #: * This checks the NUMBERS, including the ones no anchor set can reach.
        #:   ``A`` and ``T0`` enter the Kossyi form only through the combination
        #:   ``ln A - n*ln T0``, so they trade against each other EXACTLY: moving
        #:   ``T0`` from 1e4 K to 1 K and multiplying ``A`` by ``(1/1e4)^n``
        #:   leaves k identical to twelve decimals at every (Tgas, Te) point
        #:   -- measured. No anchor set of any size, at any tolerance, can
        #:   separate them. Only reading the field can.
        #:
        #: **Why this is not a second copy of the database's numbers.** The
        #: reference is the publication each entry's own shortDesc cites,
        #: transcribed once at the top of this file as ``_SVS_ARGON``,
        #: ``_BADNELL_LITHIUM`` and ``_VORONOV_LITHIUM`` -- the same constants
        #: the anchors are computed from, so there is no third transcription.
        #: Measured, every shipped entry's fields equal those published values to
        #: float noise (worst 3.7e-16 relative), so this is exact and not a
        #: tolerance exercise.
        #:
        #: **Why ``shipped_fit`` does not already do this, measured.**
        #: ``is_identical_to`` compares fields through ``ScalarQuantity.equals``,
        #: whose ``approx_equal`` accepts a 1% RELATIVE difference or a 0.01
        #: ABSOLUTE one, whichever is looser. Measured on ``BadnellRRArrhenius``:
        #: ``A`` perturbed by 0.9% compares identical, and ``B`` (0.364) by 2.0%,
        #: because the absolute term dominates a small parameter. So
        #: ``shipped_fit`` resolves a wrong ELEMENT, not parameter drift -- which
        #: is what it was added for -- and it is the anchors and this assertion,
        #: not ``shipped_fit``, that pin these numbers sharply.
        self.rate_law_parameters = dict(rate_law_parameters or {})
        #: Groups of parameters the anchors provably CANNOT separate, each an
        #: exact algebraic identity of the rate law rather than an observation
        #: that two sensitivity columns look alike.
        #:
        #: Declaring one is a statement that no anchor grid could ever do better,
        #: so it costs the entry a rank in
        #: ``test_the_anchors_separate_every_parameter_they_are_not_declared_blind_to``
        #: and is checked against the measured null direction. Anything the
        #: anchors fail to separate that is NOT declared here is a conditioning
        #: failure in the grid, and that check says so rather than absorbing it.
        self.anchor_degeneracies = tuple(tuple(group) for group in anchor_degeneracies)

    def _side(self, adjacency_lists):
        return [_structure(adjacency_list) for adjacency_list in adjacency_lists]

    def __str__(self):
        return '{0} [entry index {1!r} of library {2!r}, expected to be labelled {3!r}]'.format(
            self.name, self.entry_index, self.library, self.entry_label)

    def assert_as_shipped(self, reaction, kinetics=None, entry=None):
        """Assert that the entry selected by provenance carries what it should.

        **Every** wrong attribute is reported, not the first one. A short-circuit
        would have made this nearly useless in practice: measured on the five known
        bypasses, four of them change the kinetics class as well as the property
        that is really at issue, so a first-failure report would have said "wrong
        KINETICS CLASS" four times over and the order and response assertions would
        never have been seen to fire at all. A reader of a red run has to be able
        to see that the third-order swap is wrong in its ORDER and the Voronov
        substitution is wrong in its Te RESPONSE, because those are different
        chemistry errors with different repairs.

        A wrong value is therefore a red check that says *which* properties of the
        entry moved -- which is the thing the attribute-keyed selector could not
        do, because there a wrong value silently sent the selection somewhere else.

        `kinetics` defaults to ``reaction.kinetics``; a caller holding a raw
        ``entry.item`` (which carries ``None`` there) passes ``entry.data``.

        `entry` is the ``Entry`` the reaction was built from. It is optional only
        because a handful of checks here construct a bare reaction with no entry
        behind it; every check that selects out of a library passes it, and
        passing it is what makes the entry's own **label** an asserted attribute.
        That is the round-84 correction: the label renders the chemistry, so a
        rewritten entry has to be refused here, by name, rather than becoming a
        selection miss.
        """
        if kinetics is None:
            kinetics = reaction.kinetics
        wrong = []

        if entry is not None:
            if entry.label != self.entry_label:
                wrong.append(
                    'wrong ENTRY LABEL -- the entry selected by provenance is labelled '
                    '{0!r}, expected {1!r}. In this database an entry label IS its '
                    'reaction string: KineticsLibrary.load parses it for the '
                    'participants, the (+M) collider and the direction, and refuses any '
                    'entry whose `reversible` field disagrees with its own arrow. So a '
                    'label that has moved means the CHEMISTRY has moved, and this is the '
                    'assertion that says so instead of the selection quietly missing.'
                    .format(entry.label, self.entry_label))
            if entry.index != self.entry_index:
                wrong.append(
                    'wrong ENTRY INDEX -- the entry answers to index {0!r} but carries '
                    'index {1!r} on itself. The dict key and the entry disagree, so the '
                    'handle this suite selects on no longer identifies the entry it '
                    'names.'.format(self.entry_index, entry.index))

        if not _same_multiset(reaction.reactants, self._side(self.reactants)):
            wrong.append(
                'wrong REACTANTS -- the entry selected by provenance is not over the '
                'species this check is about.\n    got:    {0}\n    wanted: {1}'.format(
                    [str(s) for s in reaction.reactants],
                    [str(s) for s in self._side(self.reactants)]))
        if not _same_multiset(reaction.products, self._side(self.products)):
            wrong.append(
                'wrong PRODUCTS -- the entry selected by provenance is not over the '
                'species this check is about.\n    got:    {0}\n    wanted: {1}'.format(
                    [str(s) for s in reaction.products],
                    [str(s) for s in self._side(self.products)]))
        if reaction.electrons != self.electrons:
            wrong.append('wrong ELECTRONS -- net electron count is {0}, expected {1}'.format(
                reaction.electrons, self.electrons))
        if reaction.reversible != self.reversible:
            wrong.append(
                'wrong REVERSIBLE flag -- the entry is reversible={0}, expected {1}. A '
                'recombination or ionisation channel written reversible claims its own '
                'inverse is the same elementary process run backwards, which a binding '
                'project ruling forbids: the inverse of radiative recombination is '
                'photoionisation, not this entry run the other way.'.format(
                    reaction.reversible, self.reversible))
        for flag, expected in (('duplicate', self.duplicate),
                               ('allow_pdep_route', self.allow_pdep_route),
                               ('elementary_high_p', self.elementary_high_p)):
            found = getattr(reaction, flag)
            if bool(found) != bool(expected):
                wrong.append(
                    'wrong {0} flag -- the entry declares {0}={1!r}, expected {2!r}. '
                    'These three are authorable per entry and nothing else in this '
                    'suite read them, so before this assertion existed all three could '
                    'be flipped together and every check stayed green. Each changes how '
                    'RMG treats the reaction downstream: `duplicate` suppresses the '
                    'duplicate-reaction refusal, `allow_pdep_route` lets a '
                    'pressure-dependent route coexist with it, and `elementary_high_p` '
                    'declares the rate to be the high-pressure limit of an elementary '
                    'step.'.format(flag, found, expected))

        if entry is not None:
            found_amrv = bool(getattr(entry.item, 'allow_max_rate_violation', False))
            if found_amrv != bool(self.allow_max_rate_violation):
                wrong.append(
                    'wrong allow_max_rate_violation flag on the ENTRY -- '
                    'entry.item.allow_max_rate_violation={0!r}, expected {1!r}. This is '
                    'checked on the entry, not on `reaction`: '
                    '`get_library_reactions` (rmgpy/data/kinetics/library.py, the '
                    '`else:  # pdep or standard library reaction` branch) does not '
                    'forward this field into the LibraryReaction it constructs, so '
                    '`reaction.allow_max_rate_violation` is always False regardless of '
                    'what was authored -- asserting it there would be a check that can '
                    'never fail.'.format(found_amrv, self.allow_max_rate_violation))

        if kinetics is None:
            wrong.append('NO KINETICS -- the entry carries no rate coefficient at all, '
                         'so its order, its Te response and its fit cannot be read')
        else:
            if type(kinetics).__name__ != self.kinetics_class:
                wrong.append('wrong KINETICS CLASS -- the rate law is a {0}, expected a '
                             '{1}'.format(type(kinetics).__name__, self.kinetics_class))
            order = rate_order(kinetics)
            if order != self.order:
                wrong.append(
                    'wrong RATE ORDER -- the rate coefficient reads as {0}, expected {1}. '
                    'A change of order is a change of channel: second order is '
                    'A+ + e- => A, third order is A+ + 2 e- => A + e-, and they are '
                    'different reactions with different rates.'.format(
                        _describe_order(order), _describe_order(self.order)))
            response = electron_temperature_response(kinetics)
            if response != self.response:
                wrong.append(
                    'wrong Te RESPONSE -- {0}, expected {1}. A recombination rate that '
                    'rises with Te is an ionisation rate law wearing a recombination '
                    'label; electron-impact ionisation is a threshold process and '
                    'radiative recombination is not.'.format(
                        _describe_response(response), _describe_response(self.response)))
            if self.shipped_fit is not None:
                expected = self.shipped_fit()
                if not kinetics.is_identical_to(expected):
                    wrong.append(
                        'wrong FIT -- the rate law is not the one {0} gives for this '
                        'species. This is the check that sees a fit for the WRONG '
                        'ELEMENT: same class, same order, same Te response, different '
                        'atom.\n    got:    {1!r}\n    wanted: {2!r}'.format(
                            self.fit_source, kinetics, expected))
            for anchor in self.rate_anchors:
                try:
                    found = evaluate_as_the_solver_does(
                        kinetics, anchor.temperature, anchor.gas_temperature)
                except Exception as error:      # pragma: no cover - a red path
                    wrong.append(
                        'wrong RATE MAGNITUDE -- the rate law could not be evaluated at '
                        'Te = {0} K, Tgas = {1} K at all ({2}: {3}), so its magnitude is '
                        'unknown rather than right.'.format(
                            anchor.temperature, anchor.gas_temperature,
                            type(error).__name__, error))
                    continue
                if not anchor.agrees(found):
                    wrong.append(
                        'wrong RATE MAGNITUDE at Te = {0} K, Tgas = {1} K -- the entry '
                        'evaluates to {2:.6g} m^3/(mol*s), and the published fit it '
                        'claims to be gives {3:.6g} m^3/(mol*s) ({4:.3g} relative, '
                        'tolerance {5:.3g}).\n    the published form: {6}\n    this '
                        'comparison is against numbers transcribed into THIS file from '
                        'the publication, not against a second object loaded from the '
                        'same database file, so corrupting the database moves only one '
                        'side of it.'.format(
                            anchor.temperature, anchor.gas_temperature, found,
                            anchor.expected,
                            abs(found - anchor.expected) / anchor.expected,
                            anchor.rtol, anchor.source))
            # The rate law's own PARAMETERS, against the publication -- including
            # the ones no anchor can reach. See `self.rate_law_parameters` for why
            # this is not what the anchors or `shipped_fit` already do.
            for parameter, expected_value in sorted(self.rate_law_parameters.items()):
                held = getattr(kinetics, parameter, None)
                found = held if isinstance(held, float) else getattr(
                    held, 'value_si', None)
                if found is None:
                    wrong.append(
                        'wrong RATE LAW PARAMETER -- the rate law carries no {0} at '
                        'all, and the published fit gives {0} = {1!r}.'.format(
                            parameter, expected_value))
                elif not _matches_within(found, expected_value, 1e-9):
                    wrong.append(
                        'wrong RATE LAW PARAMETER -- the rate law declares {0} = {1!r} '
                        '(SI), and the publication it claims to be gives {2!r}. A '
                        'parameter can be wrong here and invisible to every rate '
                        'evaluation in this file: A and T0 enter the Kossyi form only '
                        'through ln A - n*ln T0, so they trade against each other '
                        'exactly and no anchor set can separate them.'.format(
                            parameter, found, expected_value))
            # The pressure window, which `shipped_fit` cannot see: measured,
            # `is_identical_to` ignores Pmin, Pmax and comment entirely, so a
            # reference fit built from the same shipped table compares equal to
            # one carrying a pressure range the shipped entry does not have.
            # Every entry this file names is pressure-independent -- all three
            # carry Pmin = Pmax = None, measured -- and a pressure window
            # appearing on one is a change in what the rate claims to be valid
            # over, which nothing else here would report.
            for bound in ('Pmin', 'Pmax'):
                found = getattr(kinetics, bound, None)
                if found is not None:
                    wrong.append(
                        'wrong PRESSURE WINDOW -- the rate law declares {0}={1!r}, and '
                        'every entry this file names is pressure-independent. '
                        '`is_identical_to` compares neither bound, so a shipped-fit '
                        'comparison is blind to this.'.format(bound, found))
            # The DECLARED TEMPERATURE WINDOW, which the Te-response reading cannot
            # see because it clamps itself to it -- see `self.temperature_window`.
            for bound, expected_bound in zip(('Tmin', 'Tmax'), self.temperature_window):
                found = getattr(getattr(kinetics, bound, None), 'value_si', None)
                if expected_bound is None:
                    if found is not None:
                        wrong.append(
                            'wrong TEMPERATURE WINDOW -- the rate law declares '
                            '{0}={1!r} K, and this entry is asserted to declare no '
                            '{0} at all.'.format(bound, found))
                elif not _matches_within(found, expected_bound, 1e-9):
                    wrong.append(
                        'wrong TEMPERATURE WINDOW -- the rate law declares {0}={1!r} K, '
                        'expected {2!r} K. A validity range is a claim about where this '
                        'fit may be used; the Te-response reading CLAMPS itself to it, '
                        'so moving this bound moves the points that reading is taken '
                        'at -- or is invisible to it entirely -- and shows up nowhere '
                        'else. A NON-FINITE bound is reported here too, and has to be: '
                        'measured, max(1e4, nan) is 1e4 and min(1e8, nan) is 1e8, so '
                        'the response reading clamps to the unchanged band and keeps '
                        'its answer, while production\'s own is_temperature_valid(20000) '
                        'returns False on the same object.'.format(
                            bound, found, expected_bound))

        if wrong:
            raise AssertionError('\n'.join(
                ['{0} is not as shipped -- {1} attribute(s) of the selected entry are '
                 'wrong. The entry itself is the right entry; what it CARRIES has '
                 'moved.'.format(self, len(wrong))]
                + ['  {0}'.format(item) for item in wrong]))


#: The entries this suite is about, by provenance. ``LITHIUM_CHANNEL`` maps a
#: library label to the channel of the lithium mechanism that library holds, so
#: the existing ``for label in (IONISATION, RECOMBINATION)`` loops keep their
#: shape.
#:
#: The orders, responses and flags below are what each entry is ASSERTED to
#: carry, not what it is found by. They are measured, not assumed: all three
#: shipped entries are second order (``cm^3/(molecule*s)``) -- one heavy species
#: and one electron -- and each one's Te response is read off the shipped fit.
#:
#: The ``entry_index`` values are the ``index =`` fields authored in each
#: library's ``reactions.py``, and are also the ``library.entries`` keys, measured:
#: ``PlasmaElectronImpactIonization`` holds ``{0: '[Li] => [Lip]'}`` and
#: ``PlasmaRadiativeRecombination`` holds ``{0: '[Lip] => [Li]', 1: '[Arp] =>
#: [Ar]'}``.
LITHIUM_IONISATION = ShippedEntry(
    'lithium electron-impact ionisation',
    library=IONISATION, entry_index=0, entry_label='[Li] => [Lip]',
    reactants=[LITHIUM_ATOM], products=[LITHIUM_CATION],
    order=2, response=RISES, reversible=False, electrons=1,
    kinetics_class='VoronovEIArrhenius',
    #: 1 eV to 20 keV, which is what ``voronov.yaml`` states and what
    #: ``populate_from_yaml`` converts to K with the constant restated above.
    #: Exact, not rounded: 11604.518121745585 K and 232090362.4349117 K.
    temperature_window=(1.0 / _VORONOV_KB_EV_PER_K, 2.0e4 / _VORONOV_KB_EV_PER_K),
    shipped_fit=lambda: VoronovEIArrhenius(Z=3, N=3),
    fit_source="the database's own voronov.yaml at (Z=3, N=3)",
    #: Three electron temperatures at ONE gas temperature. The Voronov form is a
    #: function of Te alone and so is its implementation, measured: k at
    #: (298.15, 2e4) and (3000, 2e4) are bit-identical. A second gas-temperature
    #: axis here would restate the first anchor and assert nothing.
    rate_anchors=[
        RateAnchor(temperature,
                   _cm3_per_molecule_to_si(_voronov_ei(temperature, **_VORONOV_LITHIUM)),
                   'Voronov (1997) ADNDT 65, 1, eq. 1 with the neutral-lithium row '
                   'A=1.39e-7 cm^3/s, dE=5.4 eV, P=0, K=0.41, X=0.438')
        for temperature in ANCHOR_TEMPERATURES],
    #: The published row itself, field by field, in SI. Measured, the shipped
    #: object carries these to 3.7e-16 relative.
    rate_law_parameters={
        'A': _cm3_per_molecule_to_si(_VORONOV_LITHIUM['A']),
        'P': _VORONOV_LITHIUM['P'],
        'X': _VORONOV_LITHIUM['X'],
        'K': _VORONOV_LITHIUM['K'],
        'dE': _VORONOV_LITHIUM['dE'],   # eV, carried on the object as a bare float
    })
LITHIUM_RECOMBINATION = ShippedEntry(
    'lithium radiative recombination',
    library=RECOMBINATION, entry_index=0, entry_label='[Lip] => [Li]',
    reactants=[LITHIUM_CATION], products=[LITHIUM_ATOM],
    order=2, response=FALLS, reversible=False, electrons=-1,
    kinetics_class='BadnellRRArrhenius',
    #: 10 K to 1e7 K, the class's own defaults, as the entry's longDesc quotes
    #: them. Wider than any condition an alkali plasma model reaches, which is
    #: why this fit -- unlike Voronov's -- is never used outside its stated range.
    temperature_window=(10.0, 1.0e7),
    shipped_fit=lambda: BadnellRRArrhenius(Z=3, N=2),
    fit_source="the database's own badnell.yaml at (Z=3, N=2)",
    #: Three electron temperatures at ONE gas temperature, for the same reason as
    #: the ionisation entry above: the Badnell form is a function of Te alone.
    rate_anchors=[
        RateAnchor(temperature,
                   _cm3_per_molecule_to_si(_badnell_rr(temperature, **_BADNELL_LITHIUM)),
                   'Badnell (2006) A&A 447, 389, eq. 1 with the Li II -> Li I row '
                   'A=8.7e-12 cm^3/s, B=0.364, T0=147.0 K, T1=7.153e6 K, C=0.1508, '
                   'T2=7.154e5 K')
        for temperature in ANCHOR_TEMPERATURES],
    #: The published row itself, field by field, in SI. Measured, the shipped
    #: object carries these to 1.8e-16 relative.
    rate_law_parameters={
        'A': _cm3_per_molecule_to_si(_BADNELL_LITHIUM['A']),
        'B': _BADNELL_LITHIUM['B'],
        'T0': _BADNELL_LITHIUM['T0'],
        'T1': _BADNELL_LITHIUM['T1'],
        'C': _BADNELL_LITHIUM['C'],
        'T2': _BADNELL_LITHIUM['T2'],
    })
#: No ``shipped_fit``: argon's recombination rate is authored inline in
#: ``reactions.py`` as a ``TwoTemperaturePlasma``, not derived from a shipped
#: table, and ``BadnellRRArrhenius(Z=18, N=17)`` does not exist -- constructing it
#: raises ``KeyError('Badnell YAML: no entry for Z=18, N=17')``, measured. A
#: reference fit here would therefore have to be a hand-copied second copy of the
#: database's numbers, and a second copy drifts.
#:
#: That left the argon rate's MAGNITUDE unchecked entirely -- ``A`` replaced by
#: ``1e100`` passed every assertion, measured -- on the entries this whole campaign
#: rests on. The anchors close it without a second copy of the database: the
#: reference is the publication the entry's own shortDesc cites, restated as a
#: formula here.
ARGON_RECOMBINATION = ShippedEntry(
    'argon radiative recombination',
    library=RECOMBINATION, entry_index=1, entry_label='[Arp] => [Ar]',
    reactants=[ARGON_CATION], products=[ARGON_ATOM],
    order=2, response=FALLS, reversible=False, electrons=-1,
    kinetics_class='TwoTemperaturePlasma',
    #: 1e4 K to 1e8 K, authored inline on the entry. Its longDesc is explicit that
    #: these are GRID MEMBERSHIP -- the span Shull & Van Steenberg themselves
    #: exercise the fit over in their Table 3 -- and NOT a quoted accuracy bound;
    #: the paper states no validity range. Asserting them here does not promote
    #: them to one. It makes moving them a decision somebody has to take on
    #: purpose, which is the only claim this assertion makes.
    temperature_window=(1.0e4, 1.0e8),
    #: SIX anchors: three electron temperatures at each of two gas temperatures.
    #: This is the only one of the three entries whose rate law is a function of
    #: the gas temperature, so it is the only one that needs -- and the only one
    #: that can use -- a second gas-temperature axis.
    #:
    #: **What the reference side is entitled to say at a second gas temperature.**
    #: Shull & Van Steenberg's power law is a function of Te only. It makes no
    #: prediction at a second Tgas on its own, and reusing a Te-only number at a
    #: second gas temperature would be asserting something the publication does
    #: not say. What is being asserted is the ENTRY's own claim, and it is a
    #: physics claim with its own justification:
    #:
    #:   The entry declares ``Ea_g = Ea_e = 0``. In the Kossyi form the gas
    #:   temperature enters only through ``exp(-Ea_g/(R*T))`` and
    #:   ``exp(Ea_e*(Te-T)/(R*T*Te))``, whose product is independent of ``T``
    #:   exactly when ``Ea_e == Ea_g``; the remaining factor is then a pure power
    #:   law in Te, which is what a radiative-recombination rate is. Radiative
    #:   recombination ``Ar+ + e- -> Ar + hv`` has no heavy-body collision partner
    #:   in it at all, so there is no gas-temperature dependence for the rate to
    #:   carry -- the only relative velocity in the problem is the electron's.
    #:   ``k`` is therefore claimed to be independent of ``Tgas``, and THAT
    #:   invariance is what these three extra anchors assert.
    #:
    #: **It does not generalise, and must not be made to.** ``PlasmaAir`` ships
    #: ``TwoTemperaturePlasma`` entries with ``Ea_e = 941.3 kJ/mol``, for which
    #: gas-temperature invariance is simply false. That is why the second axis is
    #: attached by :func:`rate_law_admits_gas_temperature_dependence` -- a measurement on the
    #: rate law in hand -- and not by a rule about the class.
    rate_anchors=[
        RateAnchor(temperature,
                   _cm3_per_molecule_to_si(_shull_van_steenberg_rr(temperature, **_SVS_ARGON)),
                   'Shull & Van Steenberg (1982) ApJS 48, 95, eq. 1 with Table 2 row '
                   'AR1, A_rad=3.77e-13 cm^3/s, eta=0.651 -- a function of Te only, '
                   'compared here at Tgas = {0} K because the entry declares '
                   'Ea_g = Ea_e = 0 and so claims k is independent of the gas '
                   'temperature'.format(gas_temperature),
                   gas_temperature=gas_temperature)
        for gas_temperature in ANCHOR_GAS_TEMPERATURES
        for temperature in ANCHOR_TEMPERATURES],
    #: The Shull & Van Steenberg row, mapped onto the Kossyi form the entry is
    #: authored in: ``A_rad`` is ``A``, ``-eta`` is ``n``, the paper's 1e4 K
    #: normalisation is ``T0``, and the two activation energies are zero because
    #: the published form has no gas-temperature dependence to carry.
    #:
    #: ``T0`` is here because it is the parameter NO anchor can reach: it enters
    #: only through ``ln A - n*ln T0``, so it trades against ``A`` exactly, and a
    #: mutation moving it from 1e4 K to 1 K with ``A`` compensated leaves every
    #: anchor identical to twelve decimals -- measured. Argon has no
    #: ``shipped_fit``, so before this round nothing in this file read it at all.
    rate_law_parameters={
        'A': _cm3_per_molecule_to_si(_SVS_ARGON['A_rad']),
        'n': -_SVS_ARGON['eta'],
        'Ea_g': 0.0,
        'Ea_e': 0.0,
        'T0': 1.0e4,
    },
    #: ``A`` and ``T0`` are one degree of freedom, not two. The Kossyi form reads
    #: ``A * (Te/T0)^n``, so in logs they appear only as ``ln A - n*ln T0``: any
    #: change to ``T0`` can be absorbed exactly by ``A``, and no anchor set at any
    #: temperatures or tolerance can tell the pair apart. Measured -- moving
    #: ``T0`` from 1e4 K to 1 K with ``A`` multiplied by ``(1/1e4)^n`` leaves k
    #: identical to twelve decimals at every anchor. This is what makes the direct
    #: parameter assertion load-bearing rather than a belt-and-braces duplicate.
    anchor_degeneracies=(('A', 'T0'),))
LITHIUM_CHANNEL = {IONISATION: LITHIUM_IONISATION, RECOMBINATION: LITHIUM_RECOMBINATION}


def _render(position, reaction, kinetics):
    """One reaction, for a failure message.

    Carries the reaction's 1-based position among the candidates searched, so
    two entries that render identically -- a duplicate-marked pair -- are still
    separately identifiable in the message, and carries the rate order and Te
    response so a reader can see at a glance what the selected entry actually
    holds.
    """
    return '#{0}  {1}  [electrons={2}, {3}, {4}, {5}]'.format(
        position, reaction, reaction.electrons,
        type(kinetics).__name__ if kinetics is not None else 'no kinetics object',
        _describe_order(rate_order(kinetics)),
        _describe_response(electron_temperature_response(kinetics)))


def _kinetics_of(reaction):
    """A reaction's own kinetics -- the default view of a candidate's rate.

    Not valid for every collection: a raw ``entry.item`` carries ``None`` here
    and holds its rate on ``entry.data`` instead (measured; see
    ``TestTheLibraryLoadLandmineIsNotReachedFromHere``).
    """
    return reaction.kinetics


def _entry_fingerprint(entry):
    """The object identities a library entry shares with every reaction built
    from it -- see :func:`select_library_reaction` for why this is the tie-back."""
    return (id(entry.data),
            tuple(sorted(id(species) for species in entry.item.reactants)),
            tuple(sorted(id(species) for species in entry.item.products)))


def _reaction_fingerprint(reaction):
    """The same identities, read off a constructed reaction."""
    return (id(reaction.kinetics),
            tuple(sorted(id(species) for species in reaction.reactants)),
            tuple(sorted(id(species) for species in reaction.products)))


def _file_digest(path):
    """sha256 of a file, so "the same library" can mean the same bytes."""
    return _file_fingerprint(path)[0]


def _file_fingerprint(path):
    """``(sha256, inode identity)`` of one file, read through ONE descriptor.

    Both halves come from the same open file -- ``handle.read()`` and
    ``os.fstat(handle.fileno())`` -- so they cannot be readings of two different
    files even if the path is replaced between them. That is the only part of this
    binding that is atomic, and it is stated plainly because the rest is not.

    **Why the sha256 alone is not enough: it cannot see an A-B-A.** A writer that
    replaces the file, lets a reader consume the replacement, and puts the original
    back leaves every later hash of that path equal to the first one. Measured, on
    this filesystem, with the content restored and the mtime forged back with
    ``os.utime``::

        before             ('9660b330', ino 103581172, ctime ...337345)
        after in-place ABA ('9660b330', ino 103581172, ctime ...519569)
            digest same: True | mtime same: True | ctime same: False
        after rename ABA   ('9660b330', ino 103581172, ctime ...627012)
            digest same: True | mtime same: True | ctime same: False
        control (no writer at all): identity unchanged

    ``st_ctime_ns`` is the discriminator: it is set to *now* by the kernel on any
    change to the inode, including a rename and including ``utime`` itself, and
    userspace cannot set it back. So it survives both shapes of the restore -- an
    in-place rewrite, and a rename-away/rename-back -- while the digest and the
    mtime survive neither honestly. ``st_dev``/``st_ino`` catch the other case, a
    different inode left at the path; ``st_size`` is belt and braces.

    **What this does NOT do.** It does not make the read atomic and it is not a
    security boundary. It is detection after the fact, like the digest, and against
    a privileged writer that can move the clock it is worth nothing.
    """
    with open(path, 'rb') as handle:
        data = handle.read()
        info = os.fstat(handle.fileno())
    return (hashlib.sha256(data).hexdigest(),
            (info.st_dev, info.st_ino, info.st_ctime_ns, info.st_size))


class _BoundLibraries:
    """The mapping :class:`LoadedLibraries` hands out in place of
    ``KineticsDatabase.libraries``.

    Deliberately a guarding mapping rather than a method, so that **every**
    existing ``database.libraries[label]`` site in this file goes through the
    identity check instead of only the two selection helpers.
    """

    def __init__(self, owner):
        self._owner = owner

    def __getitem__(self, label):
        return self._owner.library(label)

    def __iter__(self):
        return iter(self._owner.labels)

    def __len__(self):
        return len(self._owner.labels)

    def __contains__(self, label):
        return label in self._owner.labels

    def keys(self):
        return list(self._owner.labels)


class LoadedLibraries:
    """Kinetics libraries, each bound to the file it was actually loaded from.

    **The defect this answers.** A library label is not an identity.
    ``KineticsDatabase.load_libraries`` derives the label from the *directory
    name* and then does ``self.libraries[library.label] = library``
    (``rmgpy/data/kinetics/database.py:247``, and the external-path branch four
    lines above it) -- so loading a second library whose directory happens to be
    called ``PlasmaRadiativeRecombination``, from anywhere on disk, silently
    replaces the first and every later lookup by label returns the replacement.
    Two libraries that call themselves the same thing are not the same library,
    and until this class existed nothing in this suite could tell them apart:
    ``select_entry`` verified ``library.label`` and nothing else.

    **Why the library cannot answer this itself.** Measured:
    ``KineticsLibrary.load(self, path, ...)`` takes the path and does not keep it
    -- ``hasattr(library, 'path')`` is ``False`` on a loaded library, and
    ``library.name`` and ``library.label`` are both just the directory name. So
    the source path has to be captured by whoever does the loading, which is this
    class.

    **What is bound, and when.** At construction: the resolved source path of each
    label, its sha256, and the identity of the ``KineticsLibrary`` object the
    database registered for it. At every lookup: that the object registered under
    the label is still that same object, and that the source file still hashes to
    what it hashed to when it was read. The first catches a same-label library
    loaded over the top; the second catches the source being rewritten underneath
    a running check, which is precisely what the mutation drivers in
    ``evidence/`` do to a *copy* and must never do to the real tree.

    The database commit is deliberately NOT part of this. A library is identified
    by the bytes it was loaded from, which is a property of the run; pinning a
    database revision is the verifier's job, and ``evidence/verifier.sh`` does it
    there, where a legitimate database advance is a loud refusal to a human rather
    than a red test.

    **The A-B-A window, and exactly how far it is closed.** This binding is not a
    single-descriptor binding and cannot be made into one from a test file. The
    sequence is: fingerprint each source by path; call
    ``KineticsDatabase.load_libraries``, which re-opens each source **by pathname**;
    fingerprint again. A writer that replaces a source between the first
    fingerprint and the loader's open, lets the loader consume the replacement, and
    then puts the original back defeats a content hash completely -- the digest
    before, the digest after and every later digest in :meth:`library` all agree,
    over bytes that were never loaded. That race is real and this file
    *demonstrates* it rather than arguing it:
    ``test_a_source_swapped_during_the_load_does_not_pass_the_binding`` performs
    the swap at exactly the moment of the load and shows what the binding said
    before and after this repair.

    What the repair does is make the swap **detectable**: the fingerprint is
    ``(sha256, (st_dev, st_ino, st_ctime_ns, st_size))``, and ``st_ctime_ns`` moves
    on any change to the inode -- including a rename, and including ``utime``
    itself -- and cannot be set back by an unprivileged writer. Measured both ways
    in :func:`_file_fingerprint`.

    What it does NOT do, stated plainly because a mitigation described as a closure
    is worse than one described accurately:

    * it does not make the hash-load pair atomic. The loader still runs first and
      the check still runs after, so this is detection, not prevention;
    * it is not a security boundary. A writer that can move the clock, or that
      owns the filesystem, defeats it;
    * **the only true closure is an engine change and it is not made here.** It
      would be for ``KineticsLibrary.load`` to accept an already-open file object
      (or for ``KineticsDatabase.load_libraries`` to pass one down), so that the
      bytes hashed and the bytes parsed are the same descriptor. That is a change
      under ``rmgpy/``, which this branch does not touch; it is written up in
      ``docs/i134-sole-reaction-selection.md`` as an owner-gated item.
    """

    def __init__(self, libraries_path, labels):
        self.libraries_path = libraries_path
        self.labels = tuple(labels)
        assert len(set(self.labels)) == len(self.labels), (
            'the same library label was asked for twice: {0}'.format(sorted(self.labels)))
        self.sources = {}
        for label in self.labels:
            source = os.path.join(libraries_path, label, 'reactions.py')
            assert os.path.isfile(source), (
                'library {0!r} has no reactions.py at {1}, so there is nothing to bind '
                'its label to'.format(label, source))
            self.sources[label] = os.path.realpath(source)
        self._fingerprints = {label: _file_fingerprint(path)
                              for label, path in self.sources.items()}
        self._digests = {label: fingerprint[0]
                         for label, fingerprint in self._fingerprints.items()}
        self.kinetics_database = KineticsDatabase()
        self.kinetics_database.load_libraries(libraries_path, libraries=list(self.labels))
        # The A-B-A window, bracketed. `load_libraries` re-opens each source BY
        # PATHNAME, so the bytes it consumed are not the bytes hashed two lines
        # above -- a writer that swaps the file in between and puts the original
        # back afterwards leaves the digest, the mtime and the object identity all
        # agreeing over bytes that were never loaded. Re-reading the inode identity
        # after the load closes that, for the reason `_file_fingerprint` measures:
        # the restore changes st_ctime_ns even when it restores the content exactly.
        #
        # NARROWED, NOT CLOSED, and the difference matters. This is detection after
        # the fact: the loader has already run on whatever was there. A binding that
        # could not be raced at all would have to hand the loader the DESCRIPTOR
        # this class holds, and RMG's loader takes a path
        # (`KineticsDatabase.load_libraries(path, ...)` -> `KineticsLibrary.load(path)`),
        # so that is unreachable from a test file. See the class docstring for the
        # engine-side change it would take.
        self._assert_sources_did_not_move('loading them')
        missing = [label for label in self.labels
                   if label not in self.kinetics_database.libraries]
        assert not missing, (
            'asked for {0} librar(ies) and the database registered {1}; missing: '
            '{2}'.format(len(self.labels), sorted(self.kinetics_database.libraries),
                         sorted(missing)))
        self._bound = {label: self.kinetics_database.libraries[label]
                       for label in self.labels}
        self.libraries = _BoundLibraries(self)

    def _assert_sources_did_not_move(self, during, labels=None):
        """Every source file is still the file it was when it was fingerprinted.

        `during` names the interval being bracketed, so the refusal says WHEN the
        file moved rather than only that it did -- "while loading them" and "since
        it was loaded" are different failures with different repairs.
        """
        for label in (self.labels if labels is None else labels):
            path = self.sources[label]
            found = _file_fingerprint(path)
            expected = self._fingerprints[label]
            if found == expected:
                continue
            if found[0] == expected[0]:
                what = ('the CONTENT is identical and the inode is not: the file was '
                        'replaced and put back. A sha256 of the path cannot see that '
                        '-- nor can the mtime, which a writer can forge with utime -- '
                        'so the bytes this run actually read are not provably the '
                        'bytes this check names.')
            else:
                what = 'the content itself changed: {0} -> {1}'.format(expected[0], found[0])
            raise AssertionError(
                'the source of library {0!r} moved underneath this run, {1}:\n    {2}\n'
                '{3}\n  fingerprinted (dev, ino, ctime_ns, size): {4}\n'
                '  now:                                      {5}'.format(
                    label, during, path, what, expected[1], found[1]))

    def library(self, label):
        """The library bound to `label`, or an assertion naming what replaced it."""
        assert label in self._bound, (
            'library {0!r} is not one this set loaded; it loaded {1}'.format(
                label, sorted(self.labels)))
        registered = self.kinetics_database.libraries.get(label)
        if registered is not self._bound[label]:
            external = {path: name for path, name
                        in getattr(self.kinetics_database, 'external_library_labels',
                                   {}).items() if name == label}
            raise AssertionError(
                'library {0!r} is no longer the library this set loaded. It was loaded '
                'from\n    {1}\nand the object now registered under that label is a '
                'DIFFERENT object{2}.\nKineticsDatabase.load_libraries keys its registry '
                'on the directory name alone and assigns without checking '
                '(rmgpy/data/kinetics/database.py:247), so a second library whose '
                'directory is also called {0!r} -- from any checkout, any path -- '
                'replaces the first silently and every lookup by label from then on '
                'returns the replacement. A label is not an identity.'.format(
                    label, self.sources[label],
                    '' if not external else ', loaded from an external path: {0}'.format(
                        sorted(external))))
        # Reads the inode identity as well as the digest, for the reason
        # `_file_fingerprint` measures: a file replaced and put back hashes the same.
        self._assert_sources_did_not_move('since it was loaded', labels=(label,))
        assert registered.label == label, (
            'library registered under {0!r} calls itself {1!r}'.format(
                label, registered.label))
        return registered


def select_entry(library, reference):
    """The one entry of `library` that `reference` names, by provenance alone.

    :param library: a loaded :class:`KineticsLibrary`. Where the caller holds a
        :class:`LoadedLibraries`, it should be obtained through
        ``loaded.library(label)``, which is what binds the library label to the
        file it was actually loaded from -- see that class for why a label alone
        does not identify a library.
    :param reference: a :class:`ShippedEntry`; only its ``library`` and
        ``entry_index`` are read here. Neither its ``entry_label`` nor any of its
        physical attributes are consulted, on purpose: the label renders the
        chemistry in this format, so selecting on it would make a wrong-chemistry
        entry a selection miss instead of an assertion.
    :raises AssertionError: when the library is not the one named, or the index
        names no entry. The message lists every index and label the library does
        hold, which is what a human needs in order to say which entry a
        renumbered check now means.
    """
    assert library.label == reference.library, (
        'wrong library: {0} names library {1!r}, but was handed {2!r}'.format(
            reference, reference.library, library.label))
    entry = library.entries.get(reference.entry_index)
    if entry is not None:
        # ``library.entries`` is keyed on the authored index and the entry also
        # carries it, so the two can in principle disagree -- a fixture can build
        # such a library even though ``KineticsLibrary.load_entry`` cannot. The
        # disagreement is reported by ``assert_as_shipped`` as a wrong ENTRY
        # INDEX rather than raised here, so that it arrives with every other
        # wrong attribute instead of pre-empting them.
        return entry
    present = ['    index {0!r}: label {1!r}  [{2}]'.format(
        key, entry.label, type(entry.data).__name__)
        for key, entry in library.entries.items()]
    head = ('library {0!r} has no entry at index {1!r}. The index is the authored '
            '`index =` field of the entry and the key of `library.entries`; '
            'KineticsLibrary.load_entry refuses a second entry under an index already '
            'taken (library.py:680), so it is unique within a library by the loader\'s '
            'own enforcement. If the entries were RENUMBERED, this check has to be told '
            'the new number -- a renumber changes no chemistry, so nothing else in this '
            'file can work out which entry was meant. Note what has NOT happened: a '
            'rename, a reorder or a rewrite of the entry\'s chemistry all leave the '
            'index alone and are refused by assertion, by name.'.format(
                library.label, reference.entry_index))
    raise AssertionError('\n'.join(
        [head, '  wanted: {0}'.format(reference),
         '  library {0!r} holds {1} entr(ies):'.format(library.label, len(library.entries))]
        + (present or ['    (the library is empty)'])))


def select_library_reaction(library, reference, reactions=None):
    """The reaction `library` builds out of the entry `reference` names.

    Selection is by provenance, then the entry is tied back to its reaction by
    **object identity**, with no attribute comparison anywhere in the path.

    That tie-back is available because of what ``KineticsLibrary``'s loader does,
    which was probed rather than assumed (``rmgpy/data/kinetics/library.py``,
    ``get_library_reactions``): the constructed reaction *drops* ``entry.label``
    and ``entry.index`` -- so a ``LibraryReaction`` knows its library but not
    which entry it is -- while ``kinetics=entry.data`` passes the very same
    object and ``reactants=entry.item.reactants[:]`` copies the *list* but not
    the ``Species`` in it. All three branches of that function do this, including
    the ``auto_generated`` one that rewrites ``rxn.family``. Measured on the two
    libraries this suite reads: ``rxn.kinetics is entry.data`` and
    ``rxn.reactants[0] is entry.item.reactants[0]`` both hold, and repeated calls
    return fresh ``Reaction`` objects that still share those. See
    ``test_the_tie_back_is_object_identity_not_attribute_comparison``.

    :param reactions: the already-built reactions, when the caller has them;
        otherwise they are built here.
    :raises AssertionError: when the entry cannot be found (see
        :func:`select_entry`), or when the tie-back finds zero or several
        reactions -- which would mean the loader had stopped sharing objects, and
        must be loud rather than silently falling back to an attribute match.
    """
    entry = select_entry(library, reference)
    candidates = list(library.get_library_reactions() if reactions is None else reactions)
    wanted = _entry_fingerprint(entry)
    matches = [position for position, rxn in enumerate(candidates, start=1)
               if _reaction_fingerprint(rxn) == wanted]
    if len(matches) == 1:
        return candidates[matches[0] - 1]
    render = lambda position: _render(position, candidates[position - 1],
                                      _kinetics_of(candidates[position - 1]))
    lines = ['{0} reaction(s) in library {1!r} tie back to entry {2!r} by object '
             'identity, expected exactly 1. The loader shares its entry\'s kinetics '
             'object and its Species objects with every reaction it builds; if that '
             'stopped being true, this tie-back is no longer valid and must be '
             'redesigned -- it must NOT fall back to comparing '
             'attributes.'.format(len(matches), library.label, reference.entry_index),
             '  wanted: {0}'.format(reference),
             '  searched {0} reaction(s):'.format(len(candidates))]
    lines.extend('    ' + render(position) for position in range(1, len(candidates) + 1))
    raise AssertionError('\n'.join(lines))


def assert_one_reaction(reactions, where):
    """All of `reactions` are the same reaction, so ``reactions[0]`` names one.

    The other shape of the same defect: reaching a subject as ``[0]`` out of a
    collection the check does not own. It is worse than a count, because it does
    not go red when the collection grows -- it goes on silently returning *a*
    member and the check quietly changes subject.

    The repair is deliberately not ``len(reactions) == 1``, which would be the
    defect itself. What matters is not how many the producer returns but whether
    they are all one reaction.

    **What "the same reaction" means here is production's list, not a shorter
    one.** ``Reaction.is_isomorphic`` (``rmgpy/reaction.py``, the forward branch)
    conjoins four terms: the reactants, the products, ``specific_collider``, and
    the per-side electron placement counts. This guard compared three of them
    and silently dropped ``specific_collider`` -- so two pressure-dependent
    reactions differing only in their third body were "one reaction" to the
    guard and two reactions to RMG, and ``reactions[0]`` would have picked a
    collider at random while the guard certified the choice. The clause below
    closes that.

    **Both remaining terms are now production's own, not near-misses of them.**
    Until this round two of the four were paraphrases, each weaker than the thing
    it stood for, and a guard that is weaker than the predicate it licenses is a
    guard that certifies choices production would not make:

    * ``specific_collider`` was compared by graph isomorphism. Production compares
      it with ``==``, and ``Species.__eq__`` is ``self is other`` -- reference
      identity, verbatim, ``rmgpy/species.py``. So two *distinct* Species objects
      over the same molecule are one collider to isomorphism and two to RMG. The
      earlier docstring argued this direction was safe because it could only
      refuse to split reactions RMG calls one; that is backwards. It makes the
      guard *accept* -- call "one reaction" -- a pair RMG would keep as two, and
      then ``reactions[0]`` picks one of two reactions production would both have
      kept. The structural comparator that used to do this has been deleted
      rather than kept beside the operative one: two comparators for one
      question, only one of them used, is how the wrong one gets picked up by
      the next check that needs a collider compared. If a structural collider
      question is ever genuinely asked here, it should be written at that point,
      next to the reason for it.
    * the electron bookkeeping was compared as the net ``electrons`` scalar.
      Production compares :func:`get_electron_placement_counts` **per side**, and
      the whole premise of this file is that the net scalar is the comparison that
      folds ionisation onto recombination: the two shipped lithium channels have
      equal-and-opposite net counts and non-mirrored placements. A guard on the
      net scalar cannot see a pair whose placements differ but whose nets agree.
      The forward form is used -- ``counts(other) == counts(first)``, not the
      reversed pair -- because this guard asks whether the collection is one
      reaction *in one direction*, which is what licenses ``reactions[0]``.

    The consequence is stated plainly: this guard is now exactly as strict as
    production, including where production's strictness is a reference-identity
    quirk rather than a chemical statement. It will refuse a collection that holds
    two structurally identical colliders built as two objects -- and so would
    ``Reaction.is_isomorphic``, so the refusal is correct about the only question
    that matters here, which is whether ``reactions[0]`` names one thing to RMG.
    """
    assert reactions, 'nothing was generated in {0}'.format(where)
    first = reactions[0]
    first_counts = get_electron_placement_counts(first)
    for position, other in enumerate(reactions[1:], start=2):
        assert (_same_multiset(other.reactants, first.reactants)
                and _same_multiset(other.products, first.products)
                and other.specific_collider == first.specific_collider
                and get_electron_placement_counts(other) == first_counts), (
            '{0} produced {1} reactions and #{2} is not the same reaction as #1, so '
            '`reactions[0]` no longer names one thing.\n'
            '  #1:  {3} [electrons={4}, placement={5}, collider={6}]\n'
            '  #{2}:  {7} [electrons={8}, placement={9}, collider={10}]'.format(
                where, len(reactions), position, first, first.electrons, first_counts,
                _identify(first.specific_collider), other, other.electrons,
                get_electron_placement_counts(other),
                _identify(other.specific_collider)))


class _GlobalKineticsDatabase:
    """Register a kinetics database as the module-level singleton, and put back
    whatever was there before.

    ``check_for_existing_reaction`` resolves a reaction's owner through
    ``get_family_library_object``, which reads ``rmgpy.data.rmg.database`` and
    not the model's own ``kinetics_database`` attribute -- so a test driving the
    real model path has to register one. Constructing an :class:`RMGDatabase`
    does that as a side effect of ``__init__``, which means it also clobbers the
    singleton for every test that runs afterwards in the same session; hence the
    restore.
    """

    def __init__(self, kinetics_database):
        self.kinetics_database = kinetics_database
        self.saved = None

    def __enter__(self):
        self.saved = rmgpy.data.rmg.database
        database = RMGDatabase()
        database.kinetics = self.kinetics_database
        rmgpy.data.rmg.database = database
        return database

    def __exit__(self, *exc_info):
        rmgpy.data.rmg.database = self.saved
        return False


def _rate(order, response, electrons=0):
    """A minimal rate coefficient of the stated total order and Te response.

    The magnitude is meaningless and deliberately so; only the units and the
    exponent are load bearing, because those are where :func:`rate_order` and
    :func:`electron_temperature_response` read the channel from.
    ``TwoTemperaturePlasma`` is used for every combination so that the two keys
    are the *only* things that vary -- a fixture that changed the kinetics class
    at the same time would leave it open which of them the selector reacted to.

    ``response`` is required, not defaulted, for the same reason it is required on
    a :class:`ShippedEntry`: a defaulted expectation is one that is silently right
    in the fixtures and silently absent everywhere else. ``response=None`` is the
    deliberate third value -- a rate with no Te dependence at all, which
    :func:`electron_temperature_response` must call flat rather than pushing to
    one side.
    """
    units = {2: 'cm^3/(molecule*s)', 3: 'cm^6/(molecule^2*s)'}[order]
    exponent = {RISES: 1.0, FALLS: -1.0, None: 0.0}[response]
    return TwoTemperaturePlasma(A=(1.0e-13, units), n=exponent, Ea_g=(0.0, 'kJ/mol'),
                                Ea_e=(0.0, 'kJ/mol'), T0=(1e4, 'K'),
                                electrons=electrons, Tmin=(1e3, 'K'), Tmax=(1e6, 'K'))


def _reaction(reactants, products, electrons=0, family=None, kinetics=None,
              specific_collider=None, reversible=False, duplicate=False,
              allow_pdep_route=False, elementary_high_p=False,
              allow_max_rate_violation=False):
    """A reaction over shared Species objects.

    ``family`` left at ``None`` gives a plain :class:`Reaction`, which is an
    undeclared owner -- and, being a ``cdef`` class, cannot be given a ``family``
    attribute at all, so an owner has to be expressed by using a subclass that has
    one. ``LibraryReaction`` sets ``family = library``, which is how a kinetics
    library declares its electron placement on the same terms a family does.

    ``duplicate``, ``allow_pdep_route`` and ``elementary_high_p`` are here because
    they are authorable per entry and ``get_library_reactions`` copies all three
    onto the reaction it builds (``library.py:373-379``), so an assertion on the
    reaction can see them.

    ``allow_max_rate_violation`` is here for a different reason: it is accepted
    on this constructor so a fixture ``entry.item`` can author it, but
    ``get_library_reactions`` does NOT copy it onto the reaction it builds (that
    is the production defect this file names elsewhere) -- so a reaction built
    with this argument True still reads False once it comes back out of
    ``get_library_reactions``. It is here for authoring the entry side, not for
    proving anything about the reaction side.
    """
    if family is None:
        return Reaction(reactants=list(reactants), products=list(products),
                        electrons=electrons, kinetics=kinetics, reversible=reversible,
                        specific_collider=specific_collider, duplicate=duplicate,
                        allow_pdep_route=allow_pdep_route,
                        elementary_high_p=elementary_high_p,
                        allow_max_rate_violation=allow_max_rate_violation)
    return LibraryReaction(reactants=list(reactants), products=list(products),
                           electrons=electrons, library=family, reversible=reversible,
                           kinetics=kinetics, specific_collider=specific_collider,
                           duplicate=duplicate, allow_pdep_route=allow_pdep_route,
                           elementary_high_p=elementary_high_p,
                           allow_max_rate_violation=allow_max_rate_violation)


def _fixture_entry(index, label, reaction, kinetics):
    """One library entry, in the shape ``KineticsLibrary.load`` leaves them.

    ``item`` holds the reaction and ``data`` holds the rate -- which is the
    separation that makes the entry, and not the reaction, the thing a check has
    to name. ``long_desc`` is empty because a non-empty one sends
    ``get_library_reactions`` down its ``auto_generated`` branches, and this
    fixture is about the plain branch the plasma libraries take (measured:
    ``auto_generated`` is ``False`` for all five).
    """
    return Entry(index=index, label=label, item=reaction, data=kinetics, long_desc='')


def _fixture_library(label, entries):
    """A :class:`KineticsLibrary` whose ``get_library_reactions`` really runs.

    Built rather than loaded so the provenance selector is measured against the
    same production code path the database-backed checks use, without a database.
    """
    library = KineticsLibrary(label=label)
    library.auto_generated = False
    library.entries = {entry.index: entry for entry in entries}
    return library


class TestTheUndeclaredOwnerContract:
    """The contract for every reaction whose owner has never heard of electrons,
    which is all of RMG outside the plasma families and libraries.

    This is stated as a test rather than left as an inference from the suite
    passing, because it is what bounds what the repair can move. The claim is not
    "nothing broke", and it is deliberately different for the two predicates:

    * ``Reaction.is_isomorphic`` already compared the net scalar, and the
      per-side counts *reduce* to exactly that comparison for an undeclared
      owner, so its verdicts are unchanged by construction over the whole space
      of electron counts;
    * ``are_identical_species_references`` compared no electron information at
      all, so for a CHARGED reaction its verdict genuinely changes -- that is the
      repair. What bounds it is that a neutral reaction's counts are ``(0, 0)``,
      so neutral chemistry, which is all of RMG outside the charged families and
      the plasma libraries, is untouched.
    """

    ELECTRON_COUNTS = [-3, -2, -1, 0, 1, 2, 3]

    @pytest.mark.parametrize('electrons', ELECTRON_COUNTS)
    def test_counts_are_the_net_rule(self, electrons):
        rxn = _reaction([_species('[OH]')], [_species('O')], electrons=electrons)
        assert get_electron_placement_counts(rxn) == _net_rule_counts(electrons)

    def test_a_reaction_with_no_family_attribute_at_all(self):
        """A plain :class:`Reaction` -- the reactor's placement view is one -- has
        no ``family`` attribute, so the declaration lookup must not raise."""
        rxn = Reaction(reactants=[_species('[OH]')], products=[_species('O')], electrons=-1)
        assert not hasattr(rxn, 'family') or rxn.family is None
        assert get_electron_placement_counts(rxn) == (1, 0)

    def test_an_unrecognised_owner_is_undeclared(self):
        rxn = _reaction([_species('[OH]')], [_species('O')], electrons=2,
                        family='H_Abstraction')
        assert get_placement_declaration(rxn) is None
        assert get_electron_placement_counts(rxn) == (0, 2)

    @pytest.mark.parametrize('electrons1', ELECTRON_COUNTS)
    @pytest.mark.parametrize('electrons2', ELECTRON_COUNTS)
    def test_is_isomorphic_forward_verdict_is_exactly_its_old_net_comparison(
            self, electrons1, electrons2):
        """``is_isomorphic`` already compared ``self.electrons == other.electrons``
        forward. Over the whole space of pairs, the per-side comparison agrees with
        it exactly, which is what makes the change to that predicate a
        restatement for undeclared owners rather than a tightening."""
        a, b = _species('[OH]'), _species('O')
        rxn1 = _reaction([a], [b], electrons=electrons1)
        rxn2 = _reaction([a], [b], electrons=electrons2)
        assert rxn1.is_isomorphic(rxn2, either_direction=False) is (electrons1 == electrons2)

    @pytest.mark.parametrize('electrons1', ELECTRON_COUNTS)
    @pytest.mark.parametrize('electrons2', ELECTRON_COUNTS)
    def test_is_isomorphic_reverse_verdict_is_exactly_its_old_net_comparison(
            self, electrons1, electrons2):
        """Reverse identity used to be ``self.electrons == -other.electrons``.
        Same agreement, sides swapped."""
        a, b = _species('[OH]'), _species('O')
        rxn1 = _reaction([a], [b], electrons=electrons1)
        rxn2 = _reaction([b], [a], electrons=electrons2)
        assert rxn1.is_isomorphic(rxn2, either_direction=True) is (electrons1 == -electrons2)

    @pytest.mark.parametrize('electrons1', ELECTRON_COUNTS)
    @pytest.mark.parametrize('electrons2', ELECTRON_COUNTS)
    def test_the_model_predicate_is_tightened_but_only_for_charged_reactions(
            self, electrons1, electrons2):
        """``are_identical_species_references`` compared NO electron information
        before, so unlike ``is_isomorphic`` this IS a change of verdict -- for
        charged reactions. The bound on it is that two neutral reactions both have
        the counts ``(0, 0)``, so nothing about neutral chemistry moves; the
        neutral corner of this grid is the old behaviour, and the rest is the
        repair."""
        a, b = _species('[OH]'), _species('O')
        rxn1 = _reaction([a], [b], electrons=electrons1)
        rxn2 = _reaction([b], [a], electrons=electrons2)
        expected = electrons1 == -electrons2
        assert are_identical_species_references(rxn1, rxn2) is expected
        if electrons1 == 0 and electrons2 == 0:
            assert expected is True, 'two neutral reactions must still collapse'

    @pytest.mark.parametrize('owner', sorted(o for o, d in FAMILY_ELECTRON_PLACEMENT.items()
                                             if 0 in d))
    def test_one_sided_declared_owners_also_sit_on_the_net_rule(self, owner):
        """A one-sided declaration places exactly where the net rule would, so
        these owners' verdicts do not move either. Only a two-sided declaration
        can change an answer."""
        reactant_count, product_count = FAMILY_ELECTRON_PLACEMENT[owner]
        electrons = product_count - reactant_count
        rxn = _reaction([_species('[Li]')], [_species('[Li+]')],
                        electrons=electrons, family=owner)
        assert get_electron_placement_counts(rxn) == _net_rule_counts(electrons)

    def test_the_two_sided_owner_is_the_only_one_that_moves(self):
        two_sided = sorted(o for o, d in FAMILY_ELECTRON_PLACEMENT.items() if 0 not in d)
        # I-206 added the Plasma_Electron_Impact_Ionization FAMILY at (1, 2), a
        # second two-sided owner mirroring the IONISATION library's pair for the
        # same chemistry. Both "move" -- their placement differs from the net
        # rule -- which is what this test checks via the IONISATION reaction below.
        assert two_sided == [IONISATION, 'Plasma_Electron_Impact_Ionization'], two_sided
        rxn = _reaction([_species('[Li]')], [_species('[Li+]')],
                        electrons=1, family=IONISATION)
        assert get_electron_placement_counts(rxn) == (1, 2)
        assert get_electron_placement_counts(rxn) != _net_rule_counts(1)

    def test_a_reversed_declared_reaction_reports_swapped_counts(self):
        """A declaration is stated in its owner's forward orientation. A reaction
        stored the other way round carries the negated net count, and that is what
        identifies the orientation."""
        rxn = _reaction([_species('[Li+]')], [_species('[Li]')],
                        electrons=-1, family=IONISATION)
        assert get_electron_placement_counts(rxn) == (2, 1)

    def test_a_reaction_contradicting_its_owner_falls_back_rather_than_raising(self):
        """The export boundary refuses such a reaction, loudly, and should. The
        duplicate check is consulted on every comparison of every run, so it reads
        the reaction's own account of itself instead of turning a data
        inconsistency into a traceback from an unrelated place."""
        rxn = _reaction([_species('[Li]')], [_species('[Li+]')],
                        electrons=7, family=IONISATION)
        assert get_electron_placement_counts(rxn) == (0, 7)


class TestGenuineDuplicatesStillCollapse:
    """The duplicate check exists for a real reason, and the repair must not cost
    it. Two reactions that are the same reaction still have to collapse."""

    def test_identical_neutral_reactions_collapse(self):
        a, b = _species('[OH]'), _species('O')
        assert are_identical_species_references(_reaction([a], [b]), _reaction([a], [b]))

    def test_a_neutral_reaction_and_its_reverse_collapse(self):
        a, b = _species('[OH]'), _species('O')
        assert are_identical_species_references(_reaction([a], [b]), _reaction([b], [a]))

    def test_identical_charged_reactions_collapse(self):
        """Same heavy species, same owner, same placement: the same reaction."""
        a, b = _species('[Li]'), _species('[Li+]')
        rxn1 = _reaction([a], [b], electrons=1, family=IONISATION)
        rxn2 = _reaction([a], [b], electrons=1, family=IONISATION)
        assert are_identical_species_references(rxn1, rxn2)
        assert rxn1.is_isomorphic(rxn2, either_direction=True)

    def test_a_charged_reaction_and_its_true_reverse_collapse(self):
        """``A + e- => B`` and ``B => A + e-`` are one reaction seen twice, and
        RMG represents such a pair with a single object. Collapsing them is
        correct and must survive the repair."""
        a, b = _species('[Li+]'), _species('[Li]')
        forward = _reaction([a], [b], electrons=-1, family=RECOMBINATION)
        reverse = _reaction([b], [a], electrons=1, family=RECOMBINATION)
        assert get_electron_placement_counts(forward) == (1, 0)
        assert get_electron_placement_counts(reverse) == (0, 1)
        assert are_identical_species_references(forward, reverse)
        assert forward.is_isomorphic(reverse, either_direction=True)

    def test_a_different_collider_still_separates_otherwise_identical_reactions(self):
        a, b, m = _species('[OH]'), _species('O'), _species('[He]')
        rxn1 = _reaction([a], [b])
        rxn2 = _reaction([a], [b])
        rxn2.specific_collider = m
        assert not are_identical_species_references(rxn1, rxn2)

    def test_a_repeated_reaction_from_a_declared_owner_still_collapses(self):
        """A declared owner does not make two copies of one reaction distinct."""
        a, b = _species('[Li]'), _species('[Li+]')
        first = _reaction([a], [b], electrons=1, family=IONISATION)
        second = _reaction([a], [b], electrons=1, family=IONISATION)
        assert get_electron_placement_counts(first) == get_electron_placement_counts(second)
        assert are_identical_species_references(first, second)


class TestTheSelectorFailsLoudly:
    """The selector the database-backed checks reach their subject through.

    These exist because **a selector that quietly returns the first match passes
    whether or not the selection is right**, and would therefore be a worse
    defect than the count it replaced: the count at least went red. So the two
    failure directions are checked here as first-class behaviour, on their
    messages and not merely on the fact that something was raised.

    No database: the candidates are built in-process, so what is measured is the
    selector and not the shipped libraries.
    """

    def _ionisation_library(self):
        """The ionisation library's one shipped entry, rebuilt in process."""
        return _fixture_library(IONISATION, [
            _fixture_entry(0, '[Li] => [Lip]',
                           _reaction([_structure(LITHIUM_ATOM, 'Li')],
                                     [_structure(LITHIUM_CATION, 'Lip')],
                                     electrons=1, family=IONISATION),
                           VoronovEIArrhenius(Z=3, N=3)),
        ])

    def _recombination_library(self):
        """The recombination library's two shipped entries, rebuilt in process:
        lithium at index 0, argon at index 1 -- the growth that broke the old
        ``reactions[0]`` fetch."""
        return _fixture_library(RECOMBINATION, [
            _fixture_entry(0, '[Lip] => [Li]',
                           _reaction([_structure(LITHIUM_CATION, 'Lip')],
                                     [_structure(LITHIUM_ATOM, 'Li')],
                                     electrons=-1, family=RECOMBINATION),
                           BadnellRRArrhenius(Z=3, N=2)),
            _fixture_entry(1, '[Arp] => [Ar]',
                           _reaction([_structure(ARGON_CATION, 'Arp')],
                                     [_structure(ARGON_ATOM, 'Ar')],
                                     electrons=-1, family=RECOMBINATION),
                           _rate(2, FALLS, electrons=-1)),
        ])

    def test_the_named_entry_is_the_one_selected(self):
        """The control, and the reason the failures below are not vacuous."""
        library = self._recombination_library()
        found = select_library_reaction(library, LITHIUM_RECOMBINATION)
        assert found.kinetics is library.entries[0].data
        assert found.reactants[0] is library.entries[0].item.reactants[0]

    def test_the_selection_is_not_position_dependent(self):
        """Argon sits second here and first in the reversed dict; the selector
        must return the same entry either way. This is what an index-0 fetch
        cannot do, and it is also what a *reorder* of ``reactions.py`` must not
        be able to change."""
        library = self._recombination_library()
        argon = library.entries[1]
        assert select_library_reaction(library, ARGON_RECOMBINATION).kinetics is argon.data

        reordered = _fixture_library(RECOMBINATION, [library.entries[1], library.entries[0]])
        assert list(reordered.entries.values())[0] is argon, 'the fixture did not reorder'
        assert select_library_reaction(reordered, ARGON_RECOMBINATION).kinetics is argon.data

    def test_a_rename_is_still_selected_and_fails_as_an_assertion(self):
        """The round-84 correction, stated as the test that inverts.

        A rename -- which in this format is always also a chemistry change, because
        the label IS the reaction string -- used to be a **selection miss**: the
        handle was the label, so a renamed entry could not be found and no
        assertion about it ever ran. Now the handle is the index, the entry is
        still selected, and the label is refused by name along with whatever moved
        with it.
        """
        library = self._recombination_library()
        renamed = _fixture_library(RECOMBINATION, [
            _fixture_entry(0, '[Lip] + [e] => [Li]', library.entries[0].item,
                           library.entries[0].data),
            library.entries[1],
        ])
        entry = select_entry(renamed, LITHIUM_RECOMBINATION)
        assert entry is renamed.entries[0], 'the rename moved the selection'

        with pytest.raises(AssertionError) as raised:
            LITHIUM_RECOMBINATION.assert_as_shipped(
                select_library_reaction(renamed, LITHIUM_RECOMBINATION), entry=entry)
        message = str(raised.value)
        assert 'wrong ENTRY LABEL' in message
        assert "'[Lip] + [e] => [Li]'" in message, 'the message does not name what was found'
        assert "'[Lip] => [Li]'" in message, 'the message does not name what was wanted'

    def test_a_renumber_is_refused_loudly_rather_than_redirecting(self):
        """The cost of the index handle, paid where it can be seen.

        ``index =`` is the handle, so a renumber genuinely does move it. It moves
        it to *nothing*: the index is absent and the selection refuses, naming
        every index and label the library does hold. That is the one edit this
        handle cannot follow, and it is deliberately the edit that changes no
        chemistry -- a renumber is rare, a chemistry change is what the assertions
        exist for, and only one of the two can be made to fail loudly by a handle.
        """
        library = self._recombination_library()
        lithium = library.entries[0]
        renumbered = _fixture_library(RECOMBINATION, [
            _fixture_entry(41, lithium.label, lithium.item, lithium.data),
            _fixture_entry(7, library.entries[1].label, library.entries[1].item,
                           library.entries[1].data),
        ])
        assert set(renumbered.entries) == {7, 41}
        with pytest.raises(AssertionError) as raised:
            select_library_reaction(renumbered, LITHIUM_RECOMBINATION)
        message = str(raised.value)
        assert 'has no entry at index 0' in message
        assert str(LITHIUM_RECOMBINATION) in message, 'the message does not say what was wanted'
        assert 'index 41' in message and 'index 7' in message, (
            'the message does not name every index the library does hold:\n{0}'.format(message))
        assert "'[Lip] => [Li]'" in message, (
            'the message does not name the labels beside the indices, so a human cannot '
            'tell which renumbered entry was meant')

    def test_a_swap_of_two_indices_redirects_and_is_caught_by_assertion(self):
        """The residual risk of this handle, and why ``entry_label`` moved to the
        asserted half.

        Swapping two entries' indices is the only database edit that can silently
        *redirect* the selection: the index exists, so nothing refuses, and the
        entry found under it is a different reaction. It does not stay silent --
        every asserted attribute that differs between the two entries goes red by
        name, starting with the label. Four named attributes for this swap,
        measured below rather than claimed.
        """
        library = self._recombination_library()
        swapped = _fixture_library(RECOMBINATION, [
            _fixture_entry(0, library.entries[1].label, library.entries[1].item,
                           library.entries[1].data),
            _fixture_entry(1, library.entries[0].label, library.entries[0].item,
                           library.entries[0].data),
        ])
        entry = select_entry(swapped, LITHIUM_RECOMBINATION)
        assert entry.label == '[Arp] => [Ar]', 'the fixture did not swap the two indices'

        with pytest.raises(AssertionError) as raised:
            LITHIUM_RECOMBINATION.assert_as_shipped(
                select_library_reaction(swapped, LITHIUM_RECOMBINATION), entry=entry)
        message = str(raised.value)
        named = [what for what in ('wrong ENTRY LABEL', 'wrong REACTANTS', 'wrong PRODUCTS',
                                   'wrong FIT') if what in message]
        assert len(named) == 4, (
            'a redirect to a different reaction named only {0}:\n{1}'.format(named, message))

    def test_two_entries_with_one_label_no_longer_confuse_the_selection(self):
        """The collision the old handle had to refuse, now simply not a collision.

        Two entries of one library may legitimately share a label -- a
        duplicate-marked pair -- and under the label handle neither could be named.
        Under the index handle each is named by its own index, and the label
        becomes what it should be: an assertion about which chemistry the selected
        entry carries, which both of these satisfy.
        """
        library = self._recombination_library()
        twin = _fixture_entry(2, '[Lip] => [Li]', library.entries[0].item,
                              _rate(2, FALLS, electrons=-1))
        collided = _fixture_library(RECOMBINATION,
                                    [library.entries[0], library.entries[1], twin])
        labels = [entry.label for entry in collided.entries.values()]
        assert labels.count('[Lip] => [Li]') == 2, 'the fixture does not collide'

        assert select_entry(collided, LITHIUM_RECOMBINATION) is collided.entries[0]
        twin_reference = ShippedEntry(
            'the duplicate row', library=RECOMBINATION, entry_index=2,
            entry_label='[Lip] => [Li]',
            reactants=[LITHIUM_CATION], products=[LITHIUM_ATOM],
            order=2, response=FALLS, reversible=False, electrons=-1,
            kinetics_class='TwoTemperaturePlasma',
            temperature_window=(1.0e3, 1.0e6))   # `_rate`'s own window
        assert select_entry(collided, twin_reference) is twin

    def test_the_wrong_library_is_refused_rather_than_searched(self):
        """Provenance is a *pair*. Handing the recombination reference an
        ionisation library is a programming error in the check, not a search that
        happens to find nothing."""
        with pytest.raises(AssertionError) as raised:
            select_library_reaction(self._ionisation_library(), LITHIUM_RECOMBINATION)
        assert 'wrong library' in str(raised.value)

    def test_the_tie_back_is_object_identity_not_attribute_comparison(self):
        """The completeness argument, made into a check.

        Two entries of one library that agree on **every attribute anything can
        read** -- same structures, same direction, same net electrons, same
        reversibility, same kinetics class, same coefficients, therefore same
        order and same Te response -- and differ only in their labels. No
        attribute key of any width can separate these; provenance separates them
        trivially, and the reaction each entry produces is tied back to *its own*
        entry by object identity.

        This is also the guard on the premise the tie-back rests on. If
        ``KineticsLibrary.get_library_reactions`` ever stops sharing its entry's
        kinetics object and Species objects with the reaction it builds, this goes
        red -- and the repair is to redesign the tie-back, not to fall back to
        comparing attributes.
        """
        def one(index, label):
            return _fixture_entry(
                index, label,
                _reaction([_structure(LITHIUM_CATION, 'Lip')],
                          [_structure(LITHIUM_ATOM, 'Li')],
                          electrons=-1, family=RECOMBINATION),
                BadnellRRArrhenius(Z=3, N=2))

        first, second = one(0, '[Lip] => [Li]'), one(1, '[Lip] => [Li] (duplicate row)')
        library = _fixture_library(RECOMBINATION, [first, second])

        # The premise: they are indistinguishable by attribute.
        assert str(first.item) == str(second.item)
        assert first.item.electrons == second.item.electrons
        assert first.item.reversible == second.item.reversible
        assert type(first.data) is type(second.data)
        assert first.data.is_identical_to(second.data), (
            'the two fixture rates differ, so this is not the indistinguishable case')
        assert rate_order(first.data) == rate_order(second.data)
        assert (electron_temperature_response(first.data)
                == electron_temperature_response(second.data))
        assert (get_electron_placement_counts(first.item)
                == get_electron_placement_counts(second.item))

        # And the tie-back separates them anyway, by identity alone.
        reactions = library.get_library_reactions()
        assert len(reactions) == 2
        assert select_library_reaction(library, LITHIUM_RECOMBINATION).kinetics is first.data
        other = ShippedEntry('the duplicate row', library=RECOMBINATION,
                             entry_index=1,
                             entry_label='[Lip] => [Li] (duplicate row)',
                             reactants=[LITHIUM_CATION], products=[LITHIUM_ATOM],
                             order=2, response=FALLS, reversible=False, electrons=-1,
                             kinetics_class='BadnellRRArrhenius',
                             temperature_window=(10.0, 1.0e7))
        assert select_library_reaction(library, other).kinetics is second.data

        # The properties the tie-back actually relies on, asserted directly.
        for entry, rxn in zip([first, second], reactions):
            assert rxn.kinetics is entry.data, (
                'the loader no longer shares the entry kinetics object')
            assert rxn.reactants[0] is entry.item.reactants[0], (
                'the loader no longer shares the entry Species objects')
            assert rxn.products[0] is entry.item.products[0]
            assert not hasattr(rxn, 'label') or rxn.label != entry.label, (
                'the loader now carries the entry label through, so the tie-back '
                'could be stated more directly than by object identity')

    def _recombination_carrying(self, kinetics, reversible=False, label='[Lip] => [Li]',
                                duplicate=False, allow_pdep_route=False,
                                elementary_high_p=False, allow_max_rate_violation=False):
        """The lithium radiative recombination entry, with something substituted.

        This is the shape of every one of the known bypasses: the library is the
        right library, the entry is at the right index, and something about what
        that entry *carries* is wrong. Under provenance selection the entry is
        still selected -- correctly, it IS that entry -- so the assertion is what
        has to catch it, and the assertion is what names the attribute.

        Returns the pair ``(reaction, entry)``, because the entry's own label is
        an asserted attribute and ``get_library_reactions`` does not carry it
        through onto the reaction.
        """
        library = _fixture_library(RECOMBINATION, [
            _fixture_entry(0, label,
                           _reaction([_structure(LITHIUM_CATION, 'Lip')],
                                     [_structure(LITHIUM_ATOM, 'Li')],
                                     electrons=-1, family=RECOMBINATION,
                                     reversible=reversible, duplicate=duplicate,
                                     allow_pdep_route=allow_pdep_route,
                                     elementary_high_p=elementary_high_p,
                                     allow_max_rate_violation=allow_max_rate_violation),
                           kinetics)])
        return (select_library_reaction(library, LITHIUM_RECOMBINATION),
                select_entry(library, LITHIUM_RECOMBINATION))

    def test_the_shipped_entry_passes_every_assertion(self):
        """The control. Without it every red below could be the assertions
        refusing everything."""
        reaction, entry = self._recombination_carrying(BadnellRRArrhenius(Z=3, N=2))
        LITHIUM_RECOMBINATION.assert_as_shipped(reaction, entry=entry)

    def test_the_third_order_swap_is_still_selected_and_fails_on_the_rate_order(self):
        """Bypass 1. The three-body channel ``Li+ + 2 e- => Li + e-`` written in
        the radiative entry: same species, same direction, same net electrons,
        same owner and therefore the same ``(1, 0)`` placement, same falling Te
        response. The premise block below is what makes the order non-redundant --
        nothing else in this suite separates the pair.
        """
        radiative, _ = self._recombination_carrying(_rate(2, FALLS, electrons=-1))
        three_body, three_body_entry = self._recombination_carrying(
            _rate(3, FALLS, electrons=-1))

        assert str(radiative) == str(three_body), 'even their rendering is identical'
        assert radiative.electrons == three_body.electrons == -1
        assert (get_electron_placement_counts(radiative)
                == get_electron_placement_counts(three_body)), (
            'the placement declaration already separates the channels -- it does not, '
            'because it is keyed per owner and both channels sit in one library')
        assert (electron_temperature_response(three_body.kinetics)
                == electron_temperature_response(radiative.kinetics) == FALLS), (
            'three-body recombination falls with Te just as the radiative channel '
            'does -- if this ever stops being true the independence argument needs '
            're-deriving, not patching')

        with pytest.raises(AssertionError) as raised:
            LITHIUM_RECOMBINATION.assert_as_shipped(three_body, entry=three_body_entry)
        message = str(raised.value)
        assert 'wrong RATE ORDER' in message
        assert 'rate order 3' in message and 'rate order 2' in message, (
            'the message does not say both what was found and what was wanted')

    def test_an_ionisation_rate_on_the_entry_fails_on_the_te_response(self):
        """Bypass 2. The shipped lithium *ionisation* fit pasted onto the
        recombination entry: second order, so the order cannot see it. The two
        readings are independent in both directions, which is why neither can be
        dropped."""
        rising, rising_entry = self._recombination_carrying(_rate(2, RISES, electrons=-1))
        assert rate_order(rising.kinetics) == 2, (
            'the rate order already separates this, so the response adds nothing')

        with pytest.raises(AssertionError) as raised:
            LITHIUM_RECOMBINATION.assert_as_shipped(rising, entry=rising_entry)
        message = str(raised.value)
        assert 'wrong Te RESPONSE' in message
        assert 'rises' in message and 'falls' in message

    def test_a_reversible_radiative_entry_fails_on_the_reversible_flag(self):
        """Bypass 3. ``reversible = True`` on a radiative recombination claims its
        own inverse is the same elementary channel run backwards. It is not -- the
        inverse of radiative recombination is photoionisation, which needs a
        radiation field the reactor does not have. Every attribute of the rate law
        is untouched here, which is why no rate key could ever have seen it."""
        reversed_entry, reversed_entry_entry = self._recombination_carrying(
            BadnellRRArrhenius(Z=3, N=2), reversible=True)
        assert rate_order(reversed_entry.kinetics) == 2
        assert electron_temperature_response(reversed_entry.kinetics) == FALLS
        assert reversed_entry.kinetics.is_identical_to(BadnellRRArrhenius(Z=3, N=2)), (
            'the fixture changed the rate too, so this is not the reversible-only case')

        with pytest.raises(AssertionError) as raised:
            LITHIUM_RECOMBINATION.assert_as_shipped(reversed_entry,
                                                    entry=reversed_entry_entry)
        assert 'wrong REVERSIBLE flag' in str(raised.value)

    def test_a_generic_cross_section_fails_on_the_kinetics_class(self):
        """Bypass 4. An ``ElectronCollisionPlasma`` whose tabulated cross-section
        gives a second-order, Te-falling rate -- accepted by the old key *as* the
        lithium radiative recombination. The rate law is the wrong functional form
        for this channel, and the class is what says so."""
        cross_section = ElectronCollisionPlasma(
            energies=([1.0e4, 1.0e5, 1.0e6], 'J/mol'),
            sigma=([2.0e-20, 1.0e-20, 5.0e-21], 'm^2'))
        reaction, entry = self._recombination_carrying(cross_section)
        assert rate_order(reaction.kinetics) == 2, 'the fixture is not second order'

        with pytest.raises(AssertionError) as raised:
            LITHIUM_RECOMBINATION.assert_as_shipped(reaction, entry=entry)
        assert 'wrong KINETICS CLASS' in str(raised.value)
        assert 'ElectronCollisionPlasma' in str(raised.value)

    def test_a_fit_for_the_wrong_element_fails_on_the_fit(self):
        """Bypass 5, and the one no dimensional or monotonicity check can reach:
        a rate law that is right in every readable *property* and is for the wrong
        *atom*. The reference is constructed here from the database's own
        ``badnell.yaml`` at ``(Z=3, N=2)``, so nothing is hand-copied.
        """
        wrong_element = BadnellRRArrhenius(Z=6, N=5)
        assert rate_order(wrong_element) == rate_order(BadnellRRArrhenius(Z=3, N=2))
        assert (electron_temperature_response(wrong_element)
                == electron_temperature_response(BadnellRRArrhenius(Z=3, N=2))), (
            'the two fits differ in Te response, so the response already catches this')

        reaction, entry = self._recombination_carrying(wrong_element)
        with pytest.raises(AssertionError) as raised:
            LITHIUM_RECOMBINATION.assert_as_shipped(reaction, entry=entry)
        message = str(raised.value)
        assert 'wrong FIT' in message
        assert 'badnell.yaml at (Z=3, N=2)' in message
        assert 'wrong RATE MAGNITUDE' in message, (
            'a fit for the wrong element also has the wrong magnitude, and the anchor '
            'against the published Badnell row is what says so WITHOUT reading the '
            'database a second time')

    def test_a_rescaled_rate_fails_on_the_magnitude_where_no_shipped_fit_exists(self):
        """The argon defect, and the mechanism that closes it.

        ``ARGON_RECOMBINATION`` has no ``shipped_fit``: its rate is authored inline
        rather than looked up in a table the database ships, and
        ``BadnellRRArrhenius(Z=18, N=17)`` does not exist. That left the magnitude
        entirely unasserted -- the review replaced ``A`` with ``1e100`` and every
        check passed -- on the entries this campaign rests on.

        The anchors are the answer, and they are answerable *because they do not
        read the database*: the reference is the Shull & Van Steenberg (1982) power
        law restated in this file with its published Table 2 row AR1 constants. The
        two temperatures matter: the first pins the magnitude, and the second pins
        the slope, so a rate re-scaled AND re-sloped to reproduce one point still
        fails on the other.
        """
        shipped = TwoTemperaturePlasma(
            A=(3.77e-13, 'cm^3/(molecule*s)'), n=-0.651, Ea_g=(0.0, 'kJ/mol'),
            Ea_e=(0.0, 'kJ/mol'), T0=(1e4, 'K'), electrons=-1,
            Tmin=(1e4, 'K'), Tmax=(1e8, 'K'))
        assert ARGON_RECOMBINATION.shipped_fit is None, (
            'argon grew a shipped fit, so this check is no longer measuring the '
            'no-table case it was written for')

        def argon_entry(kinetics):
            library = _fixture_library(RECOMBINATION, [
                _fixture_entry(1, '[Arp] => [Ar]',
                               _reaction([_structure(ARGON_CATION, 'Arp')],
                                         [_structure(ARGON_ATOM, 'Ar')],
                                         electrons=-1, family=RECOMBINATION),
                               kinetics)])
            return (select_library_reaction(library, ARGON_RECOMBINATION),
                    select_entry(library, ARGON_RECOMBINATION))

        # The control: the shipped numbers pass.
        control_reaction, control_entry = argon_entry(shipped)
        ARGON_RECOMBINATION.assert_as_shipped(control_reaction, entry=control_entry)

        # The review's own mutation: A replaced by 1e100, everything else untouched.
        exploded = TwoTemperaturePlasma(
            A=(1e100, 'cm^3/(molecule*s)'), n=-0.651, Ea_g=(0.0, 'kJ/mol'),
            Ea_e=(0.0, 'kJ/mol'), T0=(1e4, 'K'), electrons=-1,
            Tmin=(1e4, 'K'), Tmax=(1e8, 'K'))
        assert rate_order(exploded) == rate_order(shipped) == 2
        assert electron_temperature_response(exploded) == FALLS, (
            'the rescale changed the Te response, so the response already catches it '
            'and the magnitude assertion is not what is being measured')
        reaction, entry = argon_entry(exploded)
        with pytest.raises(AssertionError) as raised:
            ARGON_RECOMBINATION.assert_as_shipped(reaction, entry=entry)
        message = str(raised.value)
        assert 'Shull & Van Steenberg' in message, (
            'the message does not name the publication that disagrees with the entry')
        for temperature in ANCHOR_TEMPERATURES:
            assert 'wrong RATE MAGNITUDE at Te = {0} K'.format(temperature) in message, (
                'the rescale is not caught at Te = {0} K:\n{1}'.format(
                    temperature, message))

        # A pure re-slope: `n` alone moved, `A` left where it was. Caught at both
        # anchors, which the compensated re-slope below deliberately is not.
        resloped_only = TwoTemperaturePlasma(
            A=(3.77e-13, 'cm^3/(molecule*s)'), n=-0.5, Ea_g=(0.0, 'kJ/mol'),
            Ea_e=(0.0, 'kJ/mol'), T0=(1e4, 'K'), electrons=-1,
            Tmin=(1e4, 'K'), Tmax=(1e8, 'K'))
        assert electron_temperature_response(resloped_only) == FALLS, (
            'the re-slope changed the Te response, so the response already catches '
            'it and the magnitude assertion is not what is being measured')
        reaction, entry = argon_entry(resloped_only)
        with pytest.raises(AssertionError) as raised:
            ARGON_RECOMBINATION.assert_as_shipped(reaction, entry=entry)
        message = str(raised.value)
        for temperature in ANCHOR_TEMPERATURES:
            assert 'wrong RATE MAGNITUDE at Te = {0} K'.format(temperature) in message, (
                'the re-slope is not caught at Te = {0} K:\n{1}'.format(
                    temperature, message))

        # A re-slope that is right at the first anchor and wrong at the second.
        resloped = TwoTemperaturePlasma(
            A=(3.77e-13 * 2.0 ** -0.651, 'cm^3/(molecule*s)'), n=0.0,
            Ea_g=(0.0, 'kJ/mol'), Ea_e=(0.0, 'kJ/mol'), T0=(1e4, 'K'), electrons=-1,
            Tmin=(1e4, 'K'), Tmax=(1e8, 'K'))
        low, high = ANCHOR_TEMPERATURES[0], ANCHOR_TEMPERATURES[-1]
        at_low = evaluate_as_the_solver_does(resloped, low)
        shipped_at_low = evaluate_as_the_solver_does(shipped, low)
        assert abs(at_low - shipped_at_low) <= 1e-9 * shipped_at_low, (
            'the re-sloped fixture does not agree with the shipped rate at the low '
            'anchor, so it does not test what a second anchor is for')
        reaction, entry = argon_entry(resloped)
        with pytest.raises(AssertionError) as raised:
            ARGON_RECOMBINATION.assert_as_shipped(reaction, entry=entry)
        assert 'wrong RATE MAGNITUDE at Te = {0} K'.format(high) in str(raised.value), (
            'only one anchor is doing any work:\n{0}'.format(raised.value))

        # The third independent mutation, and the one every check in this file was
        # blind to until this round: the electron activation energy. It is
        # invisible to `shipped_fit` (argon has none), to the rate order, to the
        # Te response (the rate still falls), and -- crucially -- it was invisible
        # to the anchors themselves while they read `get_rate_coefficient(T)`,
        # because that interface evaluates k(T, Te=T) and the Ea_e term is exactly
        # 1 there. See `evaluate_as_the_solver_does`.
        hot_electrons = TwoTemperaturePlasma(
            A=(3.77e-13, 'cm^3/(molecule*s)'), n=-0.651, Ea_g=(0.0, 'kJ/mol'),
            Ea_e=(100.0, 'kJ/mol'), T0=(1e4, 'K'), electrons=-1,
            Tmin=(1e4, 'K'), Tmax=(1e8, 'K'))
        assert rate_order(hot_electrons) == 2
        assert electron_temperature_response(hot_electrons) == FALLS, (
            'the Ea_e mutation changed the Te response, so the response already '
            'catches it and the magnitude assertion is not what is being measured')
        assert hot_electrons.is_identical_to(shipped) is False
        reaction, entry = argon_entry(hot_electrons)
        with pytest.raises(AssertionError) as raised:
            ARGON_RECOMBINATION.assert_as_shipped(reaction, entry=entry)
        message = str(raised.value)
        for temperature in ANCHOR_TEMPERATURES:
            assert 'wrong RATE MAGNITUDE at Te = {0} K'.format(temperature) in message, (
                'the Ea_e mutation is not caught at Te = {0} K:\n{1}'.format(
                    temperature, message))

    def test_each_plasma_rate_law_is_evaluated_through_the_evaluator_it_declares(self):
        """The per-class table in :func:`evaluate_as_the_solver_does`, as an assertion.

        That docstring states which of the four shipped plasma rate laws is a
        function of ``Te`` alone and which is a function of both temperatures, and
        that the one-temperature interface is the CORRECT evaluation for the first
        three and wrong only for the fourth. A table of that kind, written once and
        compared to nothing, is how this campaign's mirrors go stale. So it is
        checked here against the objects themselves.

        Three separate claims, and each one is load-bearing somewhere else in this
        file:

        1. every one of the four declares ``uses_electron_temperature``, which is
           what routes it away from the reactor's thermal branch at all;
        2. what the reactor returns is exactly what the class's own declared
           evaluator returns -- so :func:`evaluate_as_the_solver_does` is not
           re-deriving anything;
        3. for the three electron-temperature-only laws the one-temperature
           interface agrees with the declared evaluator to the last bit, and for
           ``TwoTemperaturePlasma`` it does not. That is the shape of the round-91
           defect: one class in four, not four in four.
        """
        electron_temperature_only = [
            VoronovEIArrhenius(Z=3, N=3),
            BadnellRRArrhenius(Z=3, N=2),
            ElectronCollisionPlasma(energies=([1.0e3, 1.0e4, 1.0e5, 1.0e6], 'J/mol'),
                                    sigma=([4.0e-20, 2.0e-20, 8.0e-21, 2.0e-21], 'm^2')),
        ]
        two_temperature = TwoTemperaturePlasma(
            A=(3.77e-13, 'cm^3/(molecule*s)'), n=-0.651, Ea_g=(0.0, 'kJ/mol'),
            Ea_e=(100.0, 'kJ/mol'), T0=(1e4, 'K'), electrons=-1,
            Tmin=(1e4, 'K'), Tmax=(1e8, 'K'))

        for kinetics in electron_temperature_only + [two_temperature]:
            assert getattr(kinetics, 'uses_electron_temperature', False) is True, (
                '{0} does not declare uses_electron_temperature, so the reactor would '
                'evaluate it at the gas temperature'.format(type(kinetics).__name__))

        for temperature in ANCHOR_TEMPERATURES:
            for kinetics in electron_temperature_only:
                declared = kinetics.get_rate_coefficient_electron_temp(temperature)
                assert evaluate_as_the_solver_does(kinetics, temperature) == declared, (
                    '{0} is not reached through its own declared electron-temperature '
                    'evaluator'.format(type(kinetics).__name__))
                assert kinetics.get_rate_coefficient(temperature) == declared, (
                    '{0} is an electron-temperature-only rate law whose one-temperature '
                    'interface no longer agrees with its declared evaluator, so the '
                    'claim that the old call was CORRECT for this class has '
                    'moved'.format(type(kinetics).__name__))

            declared = two_temperature.get_rate_coefficient_two_temp(
                SOLVER_TGAS, temperature)
            assert evaluate_as_the_solver_does(two_temperature, temperature) == declared, (
                'TwoTemperaturePlasma is not reached through get_rate_coefficient_two_temp '
                'at (SOLVER_TGAS, Te)')
            assert two_temperature.get_rate_coefficient(temperature) != declared, (
                'TwoTemperaturePlasma\'s one-temperature interface now agrees with its '
                'two-temperature evaluator, so the defect this round repaired no longer '
                'exists and the repair should be re-argued rather than kept on faith')

    def test_a_rate_law_that_is_not_a_function_of_te_is_refused_not_answered(self):
        """A thermal rate law asked "what do you do as Te rises" gets no number.

        Production evaluates anything that does not declare
        ``uses_electron_temperature`` at ``(Tgas, P)`` through the ordinary thermal
        branch. Reading such a rate off an electron-temperature axis is a category
        error, and the previous code committed it silently: it called
        ``get_rate_coefficient(Te)``, which for an ``Arrhenius`` is k at a GAS
        temperature of Te, and reported a direction.

        The answer is a named ``NOT_EVALUABLE``, not a number and not a pass.
        """
        thermal = Arrhenius(A=(1.0e13, 'cm^3/(mol*s)'), n=0.0, Ea=(50.0, 'kJ/mol'))
        assert getattr(thermal, 'uses_electron_temperature', False) is False
        with pytest.raises(TypeError) as raised:
            evaluate_as_the_solver_does(thermal, 2.0e4)
        assert 'uses_electron_temperature' in str(raised.value)

        response = electron_temperature_response(thermal)
        assert isinstance(response, Undetermined) and response.kind is NOT_EVALUABLE, response
        assert 'Arrhenius' in response.detail

    def test_the_one_temperature_interface_cancels_the_electron_activation_energy(self):
        """Why this file evaluates through the reactor and not through ``k(T)``.

        ``TwoTemperaturePlasma`` is the only shipped plasma rate law that is a
        function of BOTH temperatures. Its ``get_rate_coefficient(T)`` is
        documented as ``k(T, Te=T)``, and at ``T == Te`` its electron-activation
        term ``exp(Ea_e (Te - T) / (R T Te))`` is identically 1 -- so on that
        interface ``Ea_e`` has no effect whatsoever, for any value.

        This is asserted rather than described because it is the reason the Te
        response and the anchors moved to :func:`evaluate_as_the_solver_does`. If
        the one-temperature interface ever stops collapsing this way, the reason
        stops being true and this check says so.
        """
        common = dict(A=(3.77e-13, 'cm^3/(molecule*s)'), n=-0.651,
                      Ea_g=(0.0, 'kJ/mol'), T0=(1e4, 'K'), electrons=-1,
                      Tmin=(1e4, 'K'), Tmax=(1e8, 'K'))
        cold = TwoTemperaturePlasma(Ea_e=(0.0, 'kJ/mol'), **common)
        hot = TwoTemperaturePlasma(Ea_e=(100.0, 'kJ/mol'), **common)

        for temperature in ANCHOR_TEMPERATURES:
            assert cold.get_rate_coefficient(temperature) == \
                hot.get_rate_coefficient(temperature), (
                'the one-temperature interface no longer cancels Ea_e at {0} K, so the '
                'premise this file\'s evaluator change rests on has moved'.format(
                    temperature))
            through_the_reactor = (evaluate_as_the_solver_does(hot, temperature)
                                   / evaluate_as_the_solver_does(cold, temperature))
            assert through_the_reactor > 1e10, (
                'the reactor\'s own evaluator does not separate the two either '
                '(ratio {0:g} at Te = {1} K), so pointing the checks at it buys '
                'nothing'.format(through_the_reactor, temperature))

    def test_the_declared_temperature_window_is_asserted_and_not_merely_clamped_to(self):
        """A narrowed validity range, which every other check absorbs.

        :func:`electron_temperature_response` intersects :data:`TE_BAND` with the
        rate law's own ``[Tmin, Tmax]``, so narrowing the window moves the two
        points the response is read at -- and the response comes back unchanged.
        That is a check adapting to the edit it would have to catch, and it is why
        the window is now asserted by name.

        Argon is the entry where this mattered: it has no ``shipped_fit``, and
        ``is_identical_to`` -- which does compare both bounds -- is what covers the
        two lithium entries.

        **The size of the narrowing decides whether anything else sees it, and
        that is measured, not assumed.** A narrowing all the way to 9e5 K squeezes
        the response's band to [9e5, 1e6] K, where the argon power law only moves
        by a factor of 0.934 -- so THAT one does go red, as ``FLAT_IN_TE``, for a
        reason that has nothing to do with the window being wrong. A narrowing to
        1e5 K, an order of magnitude on the lower bound, leaves the response
        exactly ``FALLS`` and every other check green. The mild edit is the
        invisible one, which is the wrong way round for a check to behave.
        """
        def argon_carrying(**window):
            kinetics = TwoTemperaturePlasma(
                A=(3.77e-13, 'cm^3/(molecule*s)'), n=-0.651, Ea_g=(0.0, 'kJ/mol'),
                Ea_e=(0.0, 'kJ/mol'), T0=(1e4, 'K'), electrons=-1, **window)
            library = _fixture_library(RECOMBINATION, [
                _fixture_entry(1, '[Arp] => [Ar]',
                               _reaction([_structure(ARGON_CATION, 'Arp')],
                                         [_structure(ARGON_ATOM, 'Ar')],
                                         electrons=-1, family=RECOMBINATION),
                               kinetics)])
            return (kinetics,
                    select_library_reaction(library, ARGON_RECOMBINATION),
                    select_entry(library, ARGON_RECOMBINATION))

        # The control: the shipped window passes.
        _, control_reaction, control_entry = argon_carrying(
            Tmin=(1e4, 'K'), Tmax=(1e8, 'K'))
        ARGON_RECOMBINATION.assert_as_shipped(control_reaction, entry=control_entry)

        cases = [
            ('Tmin', '100000', dict(Tmin=(1e5, 'K'), Tmax=(1e8, 'K'))),
            # A WIDENED upper bound the response cannot see at all: TE_BAND caps
            # the sweep at 1e6 K, so moving Tmax above that changes no evaluation.
            ('Tmax', '1000000000.0', dict(Tmin=(1e4, 'K'), Tmax=(1e9, 'K'))),
        ]
        for bound, rendered, window in cases:
            kinetics, reaction, entry = argon_carrying(**window)
            assert electron_temperature_response(kinetics) == FALLS, (
                'moving {0} already shows up in the Te response, so this case is no '
                'longer measuring the hole it was written for'.format(bound))
            with pytest.raises(AssertionError) as raised:
                ARGON_RECOMBINATION.assert_as_shipped(reaction, entry=entry)
            message = str(raised.value)
            assert 'wrong TEMPERATURE WINDOW' in message
            assert bound in message and rendered in message, (
                'the message does not name the bound that moved:\n{0}'.format(message))

    @pytest.mark.parametrize('reference', [LITHIUM_IONISATION, LITHIUM_RECOMBINATION,
                                           ARGON_RECOMBINATION])
    def test_a_declared_window_that_has_a_second_source_is_tied_back_to_it(self, reference):
        """The asserted window is a second copy wherever a shipped fit exists.

        Two of the three entries here are looked up in a shipped table, and
        ``is_identical_to`` -- measured -- compares ``Tmin`` and ``Tmax``, so
        ``shipped_fit`` already covers their windows. Writing the window out again
        in this file therefore creates a mirror, and a mirror compared only to
        itself is the failure this campaign has met repeatedly: N normative copies
        agreeing with each other and none of them with the implementation.

        So the copy is tied back. Where a second source exists it must agree with
        it; where none exists -- argon, authored inline, no table to look up -- the
        value in this file is the only statement there is, which is exactly why it
        has to be here at all.
        """
        if reference.shipped_fit is None:
            assert reference is ARGON_RECOMBINATION, (
                'a second entry lost its shipped fit, so this test no longer knows '
                'which one is the untied case')
            return
        shipped = reference.shipped_fit()
        for bound, expected_bound in zip(('Tmin', 'Tmax'), reference.temperature_window):
            found = getattr(getattr(shipped, bound, None), 'value_si', None)
            assert found == expected_bound, (
                '{0}\'s asserted {1} is {2!r} K and the fit it is looked up from '
                'declares {3!r} K, so the copy in this file has drifted from its '
                'source'.format(reference.name, bound, expected_bound, found))

    @pytest.mark.parametrize('flag', ['duplicate', 'allow_pdep_route', 'elementary_high_p'])
    def test_an_authorable_flag_fails_by_name(self, flag):
        """The three entry fields nothing here used to read.

        All three are authorable per entry, all three are copied onto the reaction
        by ``get_library_reactions``, and all three could be flipped -- together,
        measured -- with every check staying green. Each changes what RMG does with
        the reaction downstream, so each is now an assertion that names itself.
        """
        reaction, entry = self._recombination_carrying(
            BadnellRRArrhenius(Z=3, N=2), **{flag: True})
        assert getattr(reaction, flag) is True, (
            'get_library_reactions no longer carries {0!r} onto the reaction, so this '
            'assertion reads something the database cannot set'.format(flag))
        with pytest.raises(AssertionError) as raised:
            LITHIUM_RECOMBINATION.assert_as_shipped(reaction, entry=entry)
        assert 'wrong {0} flag'.format(flag) in str(raised.value)

    def test_allow_max_rate_violation_is_dropped_by_get_library_reactions(self):
        """This asserts a PRODUCTION DEFECT, on purpose, so it is a notification
        and not a regression.

        Unlike ``duplicate``, ``allow_pdep_route`` and ``elementary_high_p``,
        ``allow_max_rate_violation`` is NOT copied onto the reaction that
        ``get_library_reactions`` builds. In ``rmgpy/data/kinetics/library.py``,
        inside ``get_library_reactions``, the ``else:  # pdep or standard library
        reaction`` branch constructs its ``LibraryReaction(...)`` without
        forwarding ``entry.item.allow_max_rate_violation`` -- so an entry authored
        with this flag True still comes back out of the library reading False on
        the reaction, even though the entry itself still reads True. This is a
        different ticket's fix, and the worker holding ``rmgpy/`` owns it; this
        suite is not allowed to touch that file.

        When that ticket lands, `get_library_reactions` will start forwarding the
        flag, `reaction.allow_max_rate_violation` below will become True, and THIS
        TEST WILL GO RED. That is the intended outcome: its failure is the
        notification the fix landed, and the fix should be followed by inverting
        this assertion (to ``is True``) rather than deleting the test.
        """
        reaction, entry = self._recombination_carrying(
            BadnellRRArrhenius(Z=3, N=2), allow_max_rate_violation=True)
        assert entry.item.allow_max_rate_violation is True, (
            'the fixture did not author the flag on the entry the way this test '
            'thinks it did')
        assert reaction.allow_max_rate_violation is False, (
            'get_library_reactions now forwards allow_max_rate_violation onto the '
            'reaction it builds -- the production defect named in '
            'rmgpy/data/kinetics/library.py (get_library_reactions, the '
            '`else:  # pdep or standard library reaction` branch, the '
            'LibraryReaction(...) construction) is FIXED. Invert this assertion to '
            '`is True` rather than deleting the test; this red is the notification.')
        # And the audit method catches it correctly on the entry, not on the
        # (always-False) reaction -- this is item 3(a)'s regression check firing.
        with pytest.raises(AssertionError) as raised:
            LITHIUM_RECOMBINATION.assert_as_shipped(reaction, entry=entry)
        assert 'wrong allow_max_rate_violation flag on the ENTRY' in str(raised.value)

    @pytest.mark.parametrize('bound', ['Pmin', 'Pmax'])
    def test_a_pressure_window_on_a_pressure_independent_fit_fails_by_name(self, bound):
        """The field ``shipped_fit`` is measurably blind to.

        ``is_identical_to`` compares neither ``Pmin`` nor ``Pmax`` nor the comment:
        measured, a ``BadnellRRArrhenius(Z=3, N=2)`` carrying a 1-10 bar window and
        a rewritten comment compares identical to one carrying neither. So a
        reference fit built from the shipped table cannot see a pressure window
        appearing on the entry -- and a pressure window is a change in what the
        rate claims to be valid over.
        """
        fit = BadnellRRArrhenius(Z=3, N=2)
        reference = BadnellRRArrhenius(Z=3, N=2)
        setattr(fit, bound, (1.0, 'bar'))
        assert reference.is_identical_to(fit), (
            'is_identical_to has started comparing {0}, so this check is measuring a '
            'gap that has closed and should be re-derived rather than kept'.format(bound))

        reaction, entry = self._recombination_carrying(fit)
        with pytest.raises(AssertionError) as raised:
            LITHIUM_RECOMBINATION.assert_as_shipped(reaction, entry=entry)
        message = str(raised.value)
        assert 'wrong PRESSURE WINDOW' in message
        assert bound in message
        assert 'wrong FIT' not in message, (
            'the fit comparison caught it after all, which would mean this check is '
            'redundant rather than covering a blind spot')

    def test_the_mutation_an_earlier_round_named_cannot_be_constructed(self):
        """The retraction, pinned by a check rather than by a corrected comment.

        Round 78 closed by naming ``BadnellRRArrhenius(Z=18, N=17)`` -- the shipped
        argon recombination fit put on the lithium entry -- as "the nearest
        surviving mutation and the honest next target". It does not exist: the
        database's ``badnell.yaml`` has no argon row and constructing it raises.
        A whole round's stated next step rested on a rate law nobody had tried to
        build.

        It also proposed reading the fit's own ``Z``/``N`` back off the object to
        check which element it is for. That does not work either: both classes take
        ``Z`` and ``N`` as constructor arguments, use them to look the row up, and
        keep neither. A check written on that premise would have failed on its
        first line.

        Both are asserted here so that a future round finds the refutation where it
        would look for the proposal, and not only in a report.
        """
        with pytest.raises(KeyError):
            BadnellRRArrhenius(Z=18, N=17)
        for fit in (BadnellRRArrhenius(Z=3, N=2), VoronovEIArrhenius(Z=3, N=3)):
            assert not hasattr(fit, 'Z'), (
                '{0} now keeps its Z, so the species identity CAN be read straight off '
                'the object and the fit comparison could be stated more '
                'directly'.format(type(fit).__name__))
            assert not hasattr(fit, 'N')

    def test_a_rate_with_no_te_dependence_is_named_flat_not_pushed_to_a_side(self):
        """Flat is a named third answer, not a coin toss.

        A rate coefficient independent of Te gives a ratio of exactly 1.0. Reading
        the sign of ``ratio - 1`` would have made float noise decide which
        direction it was; the decisive-factor bar makes it say so instead, and the
        assertion then reports "flat" rather than picking a side.
        """
        flat, flat_entry = self._recombination_carrying(_rate(2, None, electrons=-1))
        response = electron_temperature_response(flat.kinetics)
        assert response not in (RISES, FALLS)
        assert isinstance(response, Undetermined) and response.kind is FLAT_IN_TE, response

        with pytest.raises(AssertionError) as raised:
            LITHIUM_RECOMBINATION.assert_as_shipped(flat, entry=flat_entry)
        message = str(raised.value)
        assert 'wrong Te RESPONSE' in message
        assert 'flat' in message, 'the message does not say why the rate claimed no direction'

    def test_an_entry_with_no_kinetics_fails_on_the_absent_rate(self):
        """"No rate at all" is its own named failure, not a wildcard that passes.

        The alternative -- treating an absent rate as "unconstrained, so it
        passes" -- is the defect this campaign keeps finding, because it switches
        the check off exactly where the data is thinnest.
        """
        assert electron_temperature_response(None).kind is NO_KINETICS
        assert rate_order(None).kind is NO_KINETICS
        bare, bare_entry = self._recombination_carrying(None)
        assert bare.kinetics is None
        with pytest.raises(AssertionError) as raised:
            LITHIUM_RECOMBINATION.assert_as_shipped(bare, entry=bare_entry)
        assert 'NO KINETICS' in str(raised.value)

    def test_the_rate_order_read_here_agrees_with_the_engine_s_own_reading(self):
        """Two statements of "the order of a plasma rate", and they must agree.

        The engine has one -- ``rmgpy.electron_balance.get_plasma_rate_order`` --
        and it is a ten-entry lookup over literal unit strings, restricted to the
        four plasma kinetics classes, documented in the resolver as having "a
        silent gap, not a named failure" for anything outside it. The one here
        goes through RMG's general dimensional converter, which reads any
        spelling. Keeping both is deliberate: a selection key that simply called
        the engine's function would agree with it by construction and could not
        notice the table going stale. What must not happen is the two disagreeing
        unnoticed, which is what this pins.

        **Where the dimensional route is blind, and what it now does about it.**
        The loop below covers only the classes that carry an ``A``. Measured on
        this tree, :class:`ElectronCollisionPlasma` does not:
        ``hasattr(k, 'A')`` is ``False``, because it stores a tabulated
        cross-section. Before round 78 the dimensional route therefore reported
        "no rate-coefficient units" for it -- silently, and on a class that 9
        shipped entries across three plasma libraries, plus two plasma families,
        are written in. (That count is corrected: the claim here used to be that
        those three libraries are written WHOLLY in this class, which is false for
        two of them -- see ``TestTheShippedCrossSectionsAreRead``.) Production has an answer
        for it, and it is not a guess: ``<sigma*v>`` is bimolecular by
        construction, so the order is 2. The second block asserts that the
        dimensional route really is silent there (so the delegation is not
        redundant), that the engine really does answer (so it is not a fiction),
        and that :func:`rate_order` now returns the engine's answer.
        """
        for order in (2, 3):
            kinetics = _rate(order, FALLS)
            assert rate_order(kinetics) == order
            assert get_plasma_rate_order(kinetics) == order, (
                'the engine table and the dimensional reading disagree at order '
                '{0}'.format(order))

        cross_section = ElectronCollisionPlasma(
            energies=([0.0, 1.0e5, 1.0e6], 'J/mol'),
            sigma=([0.0, 1.0e-20, 2.0e-20], 'm^2'))
        assert not hasattr(cross_section, 'A'), (
            'ElectronCollisionPlasma grew an A-factor, so the dimensional route is '
            'no longer blind here and this delegation should be re-derived')
        assert get_plasma_rate_order(cross_section) == 2
        assert rate_order(cross_section) == 2, (
            'the rate order is blind on the class three shipped plasma libraries use')

    def test_production_already_refuses_the_swapped_channel_where_it_resolves(self):
        """The invariant is not new here -- its reach is.

        ``resolve_electron_placement`` compares the owner's declared reactant
        count against ``get_plasma_rate_order`` and refuses the view by name when
        they disagree. Measured: a ``(1, 0)`` recombination carrying a
        third-order rate is refused with ``ElectronPlacementError``, so the
        engine does catch a swapped channel -- on the path that resolves a view.

        This suite reaches library entries without ever resolving one, which is
        why the swap was invisible to it (evidence/SWAP-discriminator-absent.*:
        193 passed, 2 xfailed, exit 0, against a database whose lithium
        recombination had been turned into a different reaction). The rate order
        asserted by :meth:`ShippedEntry.assert_as_shipped` is the same invariant
        stated where this suite actually touches the data, not a second opinion
        about it.
        """
        electron = Species(molecule=[Molecule().from_adjacency_list('1 e u0 p0 c-1')],
                           label='e')

        def view_of(order):
            reaction = _reaction([_structure(LITHIUM_CATION, 'Lip')],
                                 [_structure(LITHIUM_ATOM, 'Li')],
                                 electrons=-1, family=RECOMBINATION,
                                 kinetics=_rate(order, FALLS, electrons=-1))
            species = list(reaction.reactants) + list(reaction.products) + [electron]
            return resolve_electron_placement(reaction, species)

        # The control: the shipped, second-order channel resolves.
        assert len(view_of(2).reactants) == 2

        with pytest.raises(ElectronPlacementError) as raised:
            view_of(3)
        assert 'order 3' in str(raised.value)

    def test_a_wrong_entry_is_a_selection_failure_and_a_wrong_value_is_not(self):
        """The line this round drew, asserted as behaviour rather than described.

        "The library has no lithium recombination entry" and "the lithium
        recombination entry now carries a third-order rate" are different failures
        with different repairs, and they now come out of different places: the
        first from the *selector*, naming what the library does hold; the second
        from the *assertion*, naming the attribute. Under the old attribute key
        both arrived as "nothing matches", and a third case -- a coincidence that
        made some other entry match instead -- arrived as nothing at all.

        Round 84 narrowed what can reach the selector at all: the handle is the
        entry's authored index, so "a wrong entry" now means a **renumber**, and
        every chemistry change -- including a rename, which in this format is the
        same thing -- reaches the assertions instead.
        """
        three_body = _fixture_library(RECOMBINATION, [
            _fixture_entry(0, '[Lip] => [Li]',
                           _reaction([_structure(LITHIUM_CATION, 'Lip')],
                                     [_structure(LITHIUM_ATOM, 'Li')],
                                     electrons=-1, family=RECOMBINATION),
                           _rate(3, FALLS, electrons=-1))])

        # A wrong VALUE: selected, then refused by name.
        selected = select_library_reaction(three_body, LITHIUM_RECOMBINATION)
        with pytest.raises(AssertionError) as raised:
            LITHIUM_RECOMBINATION.assert_as_shipped(selected)
        assert 'wrong RATE ORDER' in str(raised.value)

        # A wrong ENTRY: never selected at all.
        with pytest.raises(AssertionError) as raised:
            select_library_reaction(three_body, ARGON_RECOMBINATION)
        assert "has no entry at index 1" in str(raised.value)
        assert "index 0: label '[Lip] => [Li]'" in str(raised.value), (
            'the refusal does not name the entries the library does hold')

    def test_the_undetermined_outcomes_are_distinguishable_without_reading_prose(self):
        """"No rate at all", "units nothing can read", "this class states no order
        by any route" and "this composite\'s components disagree" are different
        failures with different repairs, and none of them is the integer an
        assertion asks for.

        They used to be told apart only by the first word of a diagnostic string,
        so a caller that wanted to branch on the kind had to parse English and any
        rewording changed its behaviour. The ``kind`` field is now the thing;
        ``detail`` is for humans only. This check asserts the kinds, and asserts
        that the detail is *not* what carries the distinction, by never matching on
        it.
        """

        class _Opaque:
            pass

        assert rate_order(None).kind is NO_KINETICS

        unsupported = rate_order(_Opaque())
        assert unsupported.kind is UNSUPPORTED_CLASS
        assert '_Opaque' in unsupported.detail, 'the detail does not name what it could not read'

        readable = rate_order(TwoTemperaturePlasma(
            A=(1.0, 'cm^3/(molecule*s)'), n=0.0, Ea_g=(0.0, 'kJ/mol'),
            Ea_e=(0.0, 'kJ/mol'), T0=(1e4, 'K')))
        assert readable == 2, 'the control is not a control: this one must be readable'

        class _Bogus:
            class A:
                units = 'furlongs/fortnight'
        bogus = rate_order(_Bogus())
        assert bogus.kind is UNREADABLE_UNITS, bogus
        assert 'furlongs/fortnight' in bogus.detail, 'the detail does not quote the units'

        kinds = {NO_KINETICS, UNREADABLE_UNITS, UNSUPPORTED_CLASS,
                 INCONSISTENT_COMPONENTS, EMPTY_COMPOSITE, COMPONENT_WITHOUT_UNITS,
                 PRESSURE_DEPENDENT_ORDER,
                 NO_TE_WINDOW, NOT_EVALUABLE, NOT_POSITIVE, FLAT_IN_TE}
        assert len(kinds) == 11, 'two undetermined kinds share a value, so a caller ' \
                                 'cannot tell them apart'
        assert Undetermined(NO_KINETICS, 'one wording') == Undetermined(NO_KINETICS, 'another'), (
            'equality reads the detail, so a reworded message changes what callers see')
        assert Undetermined(NO_KINETICS, 'x') != Undetermined(UNSUPPORTED_CLASS, 'x')
        assert Undetermined(UNSUPPORTED_CLASS, 'x') != 2, (
            'an undetermined order compares equal to an integer, so an == assertion '
            'against an expected order would pass on a rate that states none')
        assert Undetermined(FLAT_IN_TE, 'x') != RISES and Undetermined(FLAT_IN_TE, 'x') != FALLS

    def test_the_composite_kinetics_classes_state_an_order_after_all(self):
        """``ThirdBody``, ``Chebyshev``, ``PDepArrhenius``, ``MultiArrhenius`` and
        ``MultiPDepArrhenius`` were all reported ``unsupported`` until this round.

        That was wrong in the loudest possible way: ``unsupported`` is the outcome
        that means "no check happened here", and it was being returned for classes
        that state their dimensions perfectly well -- just not on an ``A``
        attribute. 16 shipped plasma entries are ``ThirdBody``, so this was not a
        hypothetical gap.

        A composite whose components disagree is reported as such and never
        silently resolved to the first component, and a falloff form with both
        limits is reported as pressure-dependent rather than given either limit\'s
        order.
        """
        arrhenius = Arrhenius(A=(1.0e13, 'cm^3/(mol*s)'), n=0.0, Ea=(0.0, 'kJ/mol'))
        third_order = Arrhenius(A=(1.0e13, 'cm^6/(mol^2*s)'), n=0.0, Ea=(0.0, 'kJ/mol'))

        assert rate_order(ThirdBody(arrheniusLow=third_order)) == 3
        assert rate_order(Chebyshev(coeffs=[[1.0, 0.0], [0.0, 0.0]],
                                    kunits='cm^3/(mol*s)')) == 2
        pdep = PDepArrhenius(pressures=([0.1, 1.0], 'bar'),
                             arrhenius=[arrhenius, arrhenius])
        assert rate_order(pdep) == 2
        assert rate_order(MultiArrhenius(arrhenius=[arrhenius, arrhenius])) == 2
        assert rate_order(MultiPDepArrhenius(arrhenius=[pdep, pdep])) == 2

        disagreeing = MultiArrhenius(arrhenius=[arrhenius, third_order])
        assert rate_order(disagreeing).kind is INCONSISTENT_COMPONENTS
        assert 'cm^3/(mol*s)' in rate_order(disagreeing).detail

        falloff = Lindemann(arrheniusHigh=arrhenius, arrheniusLow=third_order)
        assert rate_order(falloff).kind is PRESSURE_DEPENDENT_ORDER

    def test_a_composite_does_not_inherit_the_order_of_the_half_it_could_read(self):
        """A component whose units cannot be read is not evidence about the order.

        The defect, measured on the built module before the repair::

            component alone            -> None
            composite units            -> 'm^3/(mol*s)'
            composite rate_order       -> 2
            the valid component alone  -> 2      # indistinguishable

        The loop that walked the components skipped any component returning a bare
        ``None``, so a composite half of whose rate law states no dimensions at all
        came back with the *other* half's order, under the composite's name. Nothing
        downstream could tell that from a fully-determined composite -- which means
        the shipped-entry sweep, whose entire job is to refuse an entry that states
        no order, passed it.

        Three things are asserted, and the third is the one that matters: the
        composite is undetermined; its kind is its OWN kind and not one of the five
        that already existed (a caller that branches has to be able to say "one
        component is unreadable" apart from "the components disagree"); and the
        answer differs from the readable component's, which is the exact
        indistinguishability the measurement above showed.
        """
        readable = Arrhenius(A=(1.0e10, 'm^3/(mol*s)'), n=0.0, Ea=(0.0, 'kJ/mol'))
        unitless = Arrhenius()
        assert _rate_coefficient_units(unitless) is None, (
            'the premise moved: this component is supposed to state no units at all, '
            'and it states {0!r}'.format(_rate_coefficient_units(unitless)))
        assert rate_order(readable) == 2, 'the control is not a control'

        mixed = rate_order(MultiArrhenius(arrhenius=[readable, unitless]))
        assert isinstance(mixed, Undetermined), (
            'a composite with an unreadable component reports {0!r}, which is the '
            'readable component\'s answer wearing the composite\'s name'.format(mixed))
        assert mixed.kind is COMPONENT_WITHOUT_UNITS, mixed
        assert mixed.kind is not INCONSISTENT_COMPONENTS, (
            'an unreadable component is not a disagreement between components: the '
            'repair is to the component that states nothing, not to a conflict')
        assert mixed != rate_order(readable), (
            'the composite still answers what the readable component answers')
        assert 'MultiArrhenius' in mixed.detail and 'Arrhenius' in mixed.detail

        # Order-independent: the unreadable component first is the same answer. A
        # loop that returned on the first readable component would pass one of
        # these two and fail the other.
        other_way = rate_order(MultiArrhenius(arrhenius=[unitless, readable]))
        assert other_way.kind is COMPONENT_WITHOUT_UNITS, other_way

        # And the nested case, which is the one a real database reaches: a
        # MultiPDepArrhenius whose inner PDepArrhenius holds the unreadable part.
        nested = MultiPDepArrhenius(arrhenius=[
            PDepArrhenius(pressures=([0.1, 1.0], 'bar'), arrhenius=[readable, unitless])])
        assert rate_order(nested).kind is COMPONENT_WITHOUT_UNITS, rate_order(nested)

    def test_components_are_compared_by_order_and_not_by_how_the_units_are_spelt(self):
        """``m^3/(mol*s)`` and ``cm^3/(mol*s)`` are the same order, measured.

        The component-consistency test compared unit STRINGS, so a composite whose
        components were authored in different but dimensionally identical units was
        reported :data:`INCONSISTENT_COMPONENTS` -- the outcome that means "no check
        happened here" -- on a rate whose order is perfectly determinate. Measured
        before the repair::

            rate_order(MultiArrhenius[m^3/(mol*s), cm^3/(mol*s)])
              -> inconsistent-components: ... 2 different sets of units ...
            get_reaction_order_from_rate_coefficient_units('m^3/(mol*s)')   -> 2
            get_reaction_order_from_rate_coefficient_units('cm^3/(mol*s)')  -> 2

        Both halves are asserted here: the same order in two spellings agrees, and
        genuinely different orders still disagree -- because a repair that made
        everything agree would be the worse defect, and the two are told apart by
        nothing but this pair of cases.
        """
        metres = Arrhenius(A=(1.0e10, 'm^3/(mol*s)'), n=0.0, Ea=(0.0, 'kJ/mol'))
        centimetres = Arrhenius(A=(1.0e16, 'cm^3/(mol*s)'), n=0.0, Ea=(0.0, 'kJ/mol'))
        assert (get_reaction_order_from_rate_coefficient_units('m^3/(mol*s)')
                == get_reaction_order_from_rate_coefficient_units('cm^3/(mol*s)') == 2), (
            'the premise moved: RMG no longer reads these two spellings as one order')

        agreeing = rate_order(MultiArrhenius(arrhenius=[metres, centimetres]))
        assert agreeing == 2, (
            'two spellings of the same order are reported as {0!r}; "inconsistent" '
            'has to mean the ORDERS disagree, not that the strings do'.format(agreeing))

        third = Arrhenius(A=(1.0e16, 'cm^6/(mol^2*s)'), n=0.0, Ea=(0.0, 'kJ/mol'))
        disagreeing = rate_order(MultiArrhenius(arrhenius=[metres, third]))
        assert disagreeing.kind is INCONSISTENT_COMPONENTS, (
            'the control moved: components of genuinely different order must still '
            'be refused, or this repair made every composite agree')
        assert 'm^3/(mol*s)' in disagreeing.detail and 'cm^6/(mol^2*s)' in disagreeing.detail, (
            'the refusal no longer quotes what was actually written, so a reader of '
            'a red run cannot see which component to repair')
        assert 'order 2' in disagreeing.detail and 'order 3' in disagreeing.detail

        # A component whose spelling nothing can read is reported as unreadable
        # units and not as a disagreement: the two have different repairs.
        class _Bogus:
            class A:
                units = 'furlongs/fortnight'
        bogus = rate_order(MultiArrhenius(arrhenius=[metres, _Bogus()]))
        assert bogus.kind is UNREADABLE_UNITS, bogus
        assert 'furlongs/fortnight' in bogus.detail

    @pytest.mark.parametrize('cls', [MultiArrhenius, MultiPDepArrhenius, PDepArrhenius])
    def test_an_empty_composite_is_not_an_unsupported_class(self, cls):
        """A composite holding nothing is malformed data, not an unknown class.

        All three used to come back as :data:`UNSUPPORTED_CLASS` -- the outcome
        that means "this reader has never heard of this class, so no check
        happened". For these three that is false and misleading: the class is fully
        supported and the *entry* is broken, and the repair is to the entry.

        The three do not even agree on what empty looks like -- measured,
        ``PDepArrhenius()`` holds ``arrhenius=[]`` and the two ``Multi*`` classes
        hold ``arrhenius=None`` -- which is why :data:`COMPOSITE_CLASSES` keys on
        the class and not on the attribute. Parameterised over all three so that
        the two that arrive as ``None`` are covered and not only the one that
        happens to arrive as a list.
        """
        empty = rate_order(cls())
        assert empty.kind is EMPTY_COMPOSITE, (
            '{0}() is reported as {1!r}, not as an empty composite'.format(
                cls.__name__, empty.kind))
        assert cls.__name__ in empty.detail, 'the detail does not name the class'
        assert empty.kind is not UNSUPPORTED_CLASS

        class _Opaque:
            pass
        assert rate_order(_Opaque()).kind is UNSUPPORTED_CLASS, (
            'the control moved: a genuinely unknown class must still be unsupported')

    def test_the_one_reaction_guard_separates_on_the_specific_collider_too(self):
        """Production's identity test conjoins ``specific_collider``; so does this,
        and now by production's own comparison.

        ``Reaction.is_isomorphic`` compares reactants, products,
        ``specific_collider`` and the per-side electron counts. The guard compared
        three of the four, so a pair differing only in their third body passed it
        -- and ``reactions[0]`` then chose a collider by position while the guard
        said the choice did not matter.

        The fourth term is compared with ``==``, which for a Species is
        ``self is other``. The control below therefore passes the SAME object
        twice: two structurally identical argon atoms built as two objects are two
        reactions to production, and the block at the end asserts that this guard
        now agrees with it rather than being quietly weaker.
        """
        argon = _structure(ARGON_ATOM, 'Ar')
        helium = _species('[He]')

        def with_collider(collider):
            return _reaction([_structure(LITHIUM_ATOM, 'Li')],
                             [_structure(LITHIUM_CATION, 'Lip')],
                             electrons=1, specific_collider=collider)

        # The control: the same collider object is still one reaction -- so a green
        # result below is not the guard refusing everything.
        assert_one_reaction([with_collider(argon), with_collider(argon)], 'Fixture')
        assert_one_reaction([with_collider(None), with_collider(None)], 'Fixture')

        with pytest.raises(AssertionError) as raised:
            assert_one_reaction([with_collider(argon), with_collider(helium)], 'Fixture')
        message = str(raised.value)
        assert 'no longer names one thing' in message
        assert 'collider=' in message, 'the message does not name the colliders it split on'

        # Absent is a value, not a wildcard: a named third body and no third body
        # are two reactions, in both orders.
        with pytest.raises(AssertionError):
            assert_one_reaction([with_collider(argon), with_collider(None)], 'Fixture')
        with pytest.raises(AssertionError):
            assert_one_reaction([with_collider(None), with_collider(argon)], 'Fixture')

    def test_the_collider_comparison_is_production_s_own_reference_identity(self):
        """The alignment, checked against production and not against a
        paraphrase of production.

        Two argon atoms, isomorphic and distinct. ``Reaction.is_isomorphic`` is
        asked directly what it thinks; whatever it says, the guard must say the
        same. Written this way rather than as ``assert raises`` so that the check
        follows production if production's own comparison ever changes, instead of
        pinning today's answer as if it were the requirement.
        """
        one_argon = _structure(ARGON_ATOM, 'Ar')
        another_argon = _structure(ARGON_ATOM, 'Ar')
        assert one_argon.is_isomorphic(another_argon), 'the premise: these are the same species'
        assert one_argon is not another_argon, 'the premise: and they are two objects'

        def with_collider(collider):
            return _reaction([_structure(LITHIUM_ATOM, 'Li')],
                             [_structure(LITHIUM_CATION, 'Lip')],
                             electrons=1, specific_collider=collider)

        pair = [with_collider(one_argon), with_collider(another_argon)]
        production_says_one = pair[0].is_isomorphic(pair[1])
        assert production_says_one is False, (
            'measured premise moved: Species.__eq__ is no longer reference identity, '
            'so this whole alignment needs re-deriving')

        guard_says_one = True
        try:
            assert_one_reaction(pair, 'Fixture')
        except AssertionError:
            guard_says_one = False
        assert guard_says_one is production_says_one, (
            'the guard and Reaction.is_isomorphic disagree about whether these two '
            'are one reaction, so the guard licenses a `reactions[0]` production '
            'would not')

    def test_the_one_reaction_guard_compares_placements_not_the_net_scalar(self):
        """The second alignment: per side, as production does.

        The pair below is the premise of this whole file made into a fixture --
        the ionisation and the radiative recombination, whose net counts are equal
        and opposite and whose *placements* are ``(1, 2)`` and ``(1, 0)``, not
        mirror images. Read as one direction, a guard on the net scalar and a
        guard on the placements agree here; the case that separates them is the
        one below it, where the nets agree and the placements do not.
        """
        ionisation = _reaction([_structure(LITHIUM_ATOM, 'Li')],
                               [_structure(LITHIUM_CATION, 'Lip')],
                               electrons=1, family=IONISATION)
        undeclared = _reaction([_structure(LITHIUM_ATOM, 'Li')],
                               [_structure(LITHIUM_CATION, 'Lip')],
                               electrons=1)
        # Same species, same direction, same NET electron count -- so the old net
        # comparison called these one reaction. Their placements differ, because
        # one names a declared owner and one does not, and production splits them.
        assert ionisation.electrons == undeclared.electrons
        assert get_electron_placement_counts(ionisation) == (1, 2)
        assert get_electron_placement_counts(undeclared) == (0, 1)
        assert ionisation.is_isomorphic(undeclared) is False, (
            'production does not split these, so there is nothing to align to')

        with pytest.raises(AssertionError) as raised:
            assert_one_reaction([ionisation, undeclared], 'Fixture')
        assert 'placement=' in str(raised.value), (
            'the message does not name the placements it split on')

    def test_the_one_reaction_guard_refuses_two_different_reactions(self):
        """The other half of the repair, checked on its failure path.

        ``assert_one_reaction`` is what licenses the ``reactions[0]`` that
        reaches the attachment family's output. A guard whose failure path is
        never exercised is an assumption wearing a guard's clothes.
        """
        lithium = _reaction([_structure(LITHIUM_ATOM, 'Li')],
                            [_structure(LITHIUM_CATION, 'Lip')], electrons=1)
        argon = _reaction([_structure(ARGON_CATION, 'Arp')],
                          [_structure(ARGON_ATOM, 'Ar')], electrons=-1)
        assert_one_reaction([lithium, lithium], 'Fixture')

        with pytest.raises(AssertionError) as raised:
            assert_one_reaction([lithium, argon], 'Fixture')
        message = str(raised.value)
        assert 'no longer names one thing' in message
        assert str(lithium) in message and str(argon) in message
        assert 'Fixture produced 2 reactions' in message

    def test_the_one_reaction_guard_separates_on_the_electron_count_too(self):
        """Same heavy species, different electron bookkeeping, is a different
        reaction -- which is the premise of this entire file."""
        forward = _reaction([_structure(LITHIUM_ATOM, 'Li')],
                            [_structure(LITHIUM_CATION, 'Lip')], electrons=1)
        flipped = _reaction([_structure(LITHIUM_ATOM, 'Li')],
                            [_structure(LITHIUM_CATION, 'Lip')], electrons=-1)
        with pytest.raises(AssertionError):
            assert_one_reaction([forward, flipped], 'Fixture')

    def test_an_empty_library_says_so(self):
        with pytest.raises(AssertionError) as raised:
            select_library_reaction(_fixture_library(IONISATION, []), LITHIUM_IONISATION)
        assert 'the library is empty' in str(raised.value)

    def test_the_asserted_structures_are_adjacency_lists_and_round_trip(self):
        """The four species keys build the species they name, charge and radicals
        both. They are asserted on the selected entry, never used to find it."""
        for adjacency_list, charge, radicals in ((LITHIUM_ATOM, 0, 1), (LITHIUM_CATION, 1, 0),
                                                 (ARGON_ATOM, 0, 0), (ARGON_CATION, 1, 1)):
            molecule = _structure(adjacency_list).molecule[0]
            assert molecule.get_net_charge() == charge
            assert molecule.get_radical_count() == radicals

    @pytest.mark.xfail(strict=True, raises=AssertionError, reason=(
        'Measured on this tree: Molecule().from_smiles("[Ar+]") returns Ar(2+) -- the '
        'parser reads the radical electron as a second unit of charge -- so a '
        'SMILES-keyed structure for the argon cation is isomorphic to nothing in the '
        'shipped library. That is why the species keys in this file are adjacency lists. '
        'When this starts passing the SMILES round-trip has been repaired and this '
        'marker should go; the adjacency lists should stay either way, because they '
        'state charge and unpaired electrons separately and cannot be read two ways.'))
    def test_a_smiles_keyed_description_would_have_found_the_argon_cation(self):
        assert _species('[Ar+]').is_isomorphic(_structure(ARGON_CATION))


@pytest.mark.database
class TestChargedFamilyReactionsAreNotConfused:
    """A real family-generated pair, out of the database, with no lithium in it.

    This is the artifact that shows the defect was never a property of two plasma
    libraries meeting: it is a property of RMG's notion of reaction identity, and
    it reaches every charged reaction a family produces. The reactions here come
    from ``Plasma_Electron_Attachment`` applied to triplet dioxygen, which is what
    the family is trained on.
    """

    @classmethod
    def setup_class(cls):
        families_path = os.path.join(settings['database.directory'], 'kinetics', 'families')
        cls.database = KineticsDatabase()
        cls.database.load_recommended_families(os.path.join(families_path, 'recommended.py'))
        cls.database.load_families(families_path, families=[ATTACHMENT])
        cls.family = cls.database.families[ATTACHMENT]

    def _generated(self):
        """Every reaction the family makes from O2, regenerated per test.

        ``make_new_reaction`` replaces a reaction's lists with the model's own
        Species objects and the identity predicate compares by reference, so a
        shared reaction would carry one test's objects into the next and report a
        verdict that is an artefact of the test.
        """
        o2 = Species(molecule=[Molecule().from_adjacency_list(O2_ADJACENCY_LIST)])
        o2.generate_resonance_structures()
        reactions = self.family.generate_reactions([o2.molecule])
        assert reactions, 'the attachment family generated nothing from O2'
        return reactions

    def _attachment(self):
        """The one reaction the family makes from O2.

        Reaching a subject as ``reactions[0]`` out of a list the *family* owns is
        the same defect class as reaching it by a count over a library, with a
        worse failure mode: it does not go red when the collection grows, it goes
        on silently returning *a* reaction and the checks quietly change subject.

        Measured on this database: the family returns **two** reactions from O2,
        because the template matches each of the two equivalent oxygen atoms, and
        they are the same reaction. So the index is well defined -- and the guard
        below is what makes it well defined rather than merely lucky.
        """
        reactions = self._generated()
        assert_one_reaction(reactions, '{0} applied to O2'.format(ATTACHMENT))
        return reactions[0]

    def test_the_family_generates_exactly_one_reaction_up_to_identity(self):
        """The check that licenses the ``reactions[0]`` above.

        It is deliberately not ``len(reactions) == 1`` -- that would be the very
        defect this branch repairs, and it would be red today, because the family
        returns two. What matters is not how many the family returns but whether
        they are all one reaction; this goes red the day a genuinely different
        second recipe appears, which is the day the index stops naming a subject.
        """
        reactions = self._generated()
        assert_one_reaction(reactions, '{0} applied to O2'.format(ATTACHMENT))

    def _same_direction(self, rxn, electrons):
        return TemplateReaction(reactants=list(rxn.reactants), products=list(rxn.products),
                                electrons=electrons, family=rxn.family,
                                template=getattr(rxn, 'template', None),
                                degeneracy=rxn.degeneracy)

    def test_the_family_is_charged_and_its_reactions_carry_the_declaration(self):
        attach = self._attachment()
        assert self.family.electrons == -1
        assert attach.electrons == -1
        assert get_electron_placement_counts(attach) == (1, 0)

    def test_its_true_reverse_still_collapses(self):
        """The control. A reaction and its reverse are one reaction, and the
        check must keep saying so."""
        attach = self._attachment()
        reverse = TemplateReaction(reactants=list(attach.products),
                                   products=list(attach.reactants),
                                   electrons=-attach.electrons, family=attach.family,
                                   template=getattr(attach, 'template', None),
                                   degeneracy=attach.degeneracy)
        assert are_identical_species_references(attach, reverse)

    def test_the_same_direction_charge_reversal_survives(self):
        """``O2 + e- => O2-`` and ``O2 => O2- + e-`` are different reactions --
        different molecularity, opposite charge transfer -- and used to collapse.

        ``flipped`` is CONSTRUCTED, not generated. Measured (probe against this
        database, loading only the attachment family plus the two shipped plasma
        libraries ``PlasmaElectronImpactIonization`` /
        ``PlasmaRadiativeRecombination``): ``family.generate_reactions([o2])``
        returns 2 reactions from O2, both carrying ``electrons == -1`` (the
        family's own ``self.electrons``, fixed at family scope and stamped onto
        every reaction it produces regardless of direction --
        ``rmgpy/data/kinetics/family.py:1781``); no reaction with
        ``electrons == +1`` comes out. The two loaded libraries do not mention O2
        at all -- their entries are Li/Li+ and Ar/Ar+ only. So on this database
        neither the exercised family path nor the exercised library path can
        produce the same-direction, opposite-charge O2 state this test reads;
        the only way to get that object is ``_same_direction(attach,
        electrons=-attach.electrons)`` building it by hand.

        What this therefore proves: the identity/placement-counting machinery
        (``get_electron_placement_counts``, ``are_identical_species_references``)
        treats the two placements as distinct when handed both, as an object-level
        property. It does NOT prove the attachment family or these libraries can
        ever produce both placements for the same heavy species -- that path is
        not reachable here, and this test cannot be read as evidence that it is.
        """
        attach = self._attachment()
        flipped = self._same_direction(attach, electrons=-attach.electrons)
        assert get_electron_placement_counts(attach) == (1, 0)
        assert get_electron_placement_counts(flipped) == (0, 1)
        assert not are_identical_species_references(attach, flipped)

    def test_the_neutral_transformation_over_the_same_heavy_species_survives(self):
        """The shape that needs no plasma library at all: a neutral reaction whose
        heavy species happen to coincide with a charged one's."""
        attach = self._attachment()
        neutral = self._same_direction(attach, electrons=0)
        assert get_electron_placement_counts(neutral) == (0, 0)
        assert not are_identical_species_references(attach, neutral)

    def test_the_model_keeps_both_when_one_side_is_constructed(self):
        """``CoreEdgeReactionModel.make_new_reaction`` keeps a genuinely
        family-generated attachment AND its hand-built same-direction,
        opposite-charge counterpart as two distinct new reactions.

        Renamed from ``..._through_the_real_path``: only ``first_new`` (the
        attachment itself, via :meth:`_attachment`) comes off a real generated
        path. ``flipped`` is built by :meth:`_same_direction` with
        ``electrons=1`` -- see the measurement in
        ``test_the_same_direction_charge_reversal_survives`` above: nothing this
        family or these libraries generate on this database ever carries
        ``electrons == +1`` for O2, so ``flipped`` is not reachable through the
        model either. This test proves the model's duplicate-collapse machinery
        does not fold the two placements together when both are handed to it --
        it does not prove the model can ever be handed ``flipped`` by a real run.
        """
        with _GlobalKineticsDatabase(self.database):
            model = CoreEdgeReactionModel()
            model.kinetics_database = self.database
            _, first_new = model.make_new_reaction(self._attachment(), generate_thermo=False,
                                                   generate_kinetics=False)
            assert first_new is True
            flipped = self._same_direction(self._attachment(), electrons=1)
            _, second_new = model.make_new_reaction(flipped, generate_thermo=False,
                                                    generate_kinetics=False)
        assert second_new is True, ('the same-direction charge reversal was discarded as a '
                                    'duplicate of the attachment')


class _PlasmaLibraryFixture:
    """The two shipped plasma libraries, and selection out of them.

    Shared by the two classes below, which had identical ``setup_class`` and
    ``_database`` pairs and two divergent copies of the count-based fetch -- one
    with a message, one bare. One selector, used by both, is the point: three
    copies of a fetch is three places for the next one to rot.

    Not collected by pytest: ``pytest.ini`` collects ``*Test`` and ``Test*``
    classes, and this is neither.
    """

    @classmethod
    def setup_class(cls):
        cls.libraries_path = os.path.join(settings['database.directory'],
                                          'kinetics', 'libraries')

    def _database(self):
        """The two libraries, each bound to the file it came from.

        A :class:`LoadedLibraries`, not a bare ``KineticsDatabase``: its
        ``.libraries`` mapping guards every lookup against a same-label library
        having been loaded over the top, which a bare registry cannot do. The real
        ``KineticsDatabase`` is still available as ``.kinetics_database`` for the
        checks that hand it to ``CoreEdgeReactionModel``.
        """
        return LoadedLibraries(self.libraries_path, [IONISATION, RECOMBINATION])

    def _select(self, database, reference):
        """The :class:`LibraryReaction` built from the entry `reference` names.

        ``get_library_reactions`` is the wrapping step that attaches the owner,
        so this is the form every consumer downstream of the loader sees -- and it
        is the step that drops ``entry.label``, which is why the tie-back inside
        :func:`select_library_reaction` is by object identity.
        """
        return select_library_reaction(database.libraries[reference.library], reference)

    def _select_entry(self, database, reference):
        """The library *entry* `reference` names.

        One step earlier than :meth:`_select`: the entry still holds its reaction
        on ``item`` and its rate on ``data``, and ``entry.item.kinetics`` is
        ``None`` (measured in ``test_a_loaded_entry_carries_no_kinetics_either``).
        That is the form ``check_for_duplicates`` sees, which is what
        :class:`TestTheLibraryLoadLandmineIsNotReachedFromHere` is about.

        Under provenance this step needs no special handling at all. The old
        attribute-keyed selector had to be handed a ``kinetics_for`` override here,
        because it read the rate off the candidate and the candidate had none;
        provenance reads only the library and the entry's index, so the
        one-step-earlier incompleteness that forced the override is simply not on
        the selection path any more.
        """
        return select_entry(database.libraries[reference.library], reference)


@pytest.mark.database
class TestLithiumChargeNetworkReachesTheModel(_PlasmaLibraryFixture):
    """The real mechanism, through the real model-building path.

    ``add_reaction_library_to_edge`` is the step that dropped the sink, and
    ``CoreEdgeReactionModel`` is the model the assertion is about: the requirement
    is that the recombination is *in the model*, not merely that it loads.

    The reactions are reached by provenance rather than by index, because
    ``PlasmaRadiativeRecombination`` also carries an argon entry that has nothing
    to do with anything here, and will carry more.
    """

    def test_the_selection_discriminates_within_one_library(self):
        """The check that the selection is selecting.

        ``PlasmaRadiativeRecombination`` holds the lithium and the argon channel.
        Naming each entry returns a *different* reaction -- which is what a
        first-match or an index-0 fetch could not do, and is therefore the check
        that would go red if the selector quietly degraded into one.
        """
        database = self._database()
        lithium = self._select(database, LITHIUM_RECOMBINATION)
        argon = self._select(database, ARGON_RECOMBINATION)
        assert lithium is not argon
        assert lithium.reactants[0].is_isomorphic(_structure(LITHIUM_CATION))
        assert argon.reactants[0].is_isomorphic(_structure(ARGON_CATION))
        assert not lithium.reactants[0].is_isomorphic(_structure(ARGON_CATION))

    def test_every_named_entry_carries_what_it_is_supposed_to(self):
        """The assertion half, run end to end against the shipped database.

        This is where the rate-order and Te-response work now lives. Each entry is
        found by provenance and then checked, attribute by attribute, against what
        its channel implies -- order, Te response, reversibility, net electrons,
        kinetics class, and for the two lithium channels the shipped fit itself,
        reconstructed from the database's own ``badnell.yaml`` / ``voronov.yaml``
        by ``(Z, N)``.

        Every one of the known bypasses lands here, and each lands as a named
        attribute rather than as a selection that quietly went somewhere else.

        The entry is passed alongside the reaction because the entry's own
        **label** is one of the asserted attributes: ``get_library_reactions``
        drops it, and it is the label that carries a rewritten reaction string --
        the reversibility bypass among them -- into a named assertion instead of a
        selection miss.
        """
        database = self._database()
        for reference in (LITHIUM_IONISATION, LITHIUM_RECOMBINATION, ARGON_RECOMBINATION):
            reference.assert_as_shipped(self._select(database, reference),
                                        entry=self._select_entry(database, reference))

    def test_no_shipped_plasma_entry_label_collides(self):
        """Entry labels, measured on every plasma library the database ships.

        A label collision no longer breaks the handle -- the handle is the index,
        and the loader refuses a file that reuses one. It is kept because it is
        still a fact worth knowing about the database, and because the label is now
        an *asserted* attribute: two entries rendering the same reaction string in
        one library means two entries claiming to be the same chemistry, which is
        what ``check_for_duplicates`` exists to arbitrate.
        """
        libraries = list(self.SHIPPED_CHANNELS) + list(self.LIBRARIES_NOT_SWEPT)
        database = LoadedLibraries(self.libraries_path, libraries)
        assert set(database.libraries) == set(libraries), (
            'not every plasma library loaded, so this sweep is over fewer than it claims')
        for label in libraries:
            entries = list(database.libraries[label].entries.values())
            assert entries, '{0} loaded no entries, so this check is vacuous on it'.format(label)
            labels = [entry.label for entry in entries]
            duplicated = sorted({name for name in labels if labels.count(name) > 1})
            assert not duplicated, (
                '{0} carries {1} entr(ies) under duplicated label(s) {2}, so two entries '
                'of one library render the same chemistry'.format(
                    label, len(entries), duplicated))

    def test_a_second_library_of_the_same_name_does_not_pass_as_the_first(self, tmp_path):
        """Two libraries calling themselves the same thing are not the same library.

        Built, not argued: a byte-for-byte copy of ``PlasmaRadiativeRecombination``
        is placed in a scratch directory of the same name and loaded into the same
        ``KineticsDatabase`` by path. ``load_libraries`` takes the label from the
        directory name and assigns into ``self.libraries`` without checking
        (``rmgpy/data/kinetics/database.py:247``), so the registry entry is
        replaced.

        Two things are asserted, and the first is the finding:

        1. The overwrite really happens, and the old label check cannot see it --
           the impostor's ``label`` is identical, so ``library.label ==
           reference.library`` is as true of the replacement as of the original.
        2. :class:`LoadedLibraries` refuses it, naming the file the library it
           bound was actually loaded from.

        The copy is *identical* on purpose. An impostor with different chemistry
        would be caught by the assertions further down and would prove nothing
        about identity; this one differs from the original in nothing but where it
        came from, which is exactly the case a label cannot distinguish.
        """
        database = self._database()
        original = database.library(RECOMBINATION)
        assert original.label == RECOMBINATION

        impostor_dir = tmp_path / 'elsewhere' / RECOMBINATION
        shutil.copytree(os.path.join(self.libraries_path, RECOMBINATION), str(impostor_dir))
        database.kinetics_database.load_libraries(str(tmp_path / 'elsewhere'),
                                                  libraries=[str(impostor_dir)])
        impostor = database.kinetics_database.libraries[RECOMBINATION]

        assert impostor is not original, (
            'the second load did not replace the first, so this check is not measuring '
            'the overwrite it claims to measure')
        assert impostor.label == original.label == RECOMBINATION, (
            'the impostor does not call itself the same thing, so the label check would '
            'have caught it and nothing here is at issue')

        with pytest.raises(AssertionError) as raised:
            database.library(RECOMBINATION)
        message = str(raised.value)
        assert 'no longer the library this set loaded' in message
        assert database.sources[RECOMBINATION] in message, (
            'the refusal does not name the file the bound library came from, which is '
            'the one thing a human needs in order to tell the two apart')

        with pytest.raises(AssertionError):
            select_entry(database.libraries[RECOMBINATION], LITHIUM_RECOMBINATION)

    def test_a_source_swapped_during_the_load_does_not_pass_the_binding(
            self, tmp_path, monkeypatch):
        """The A-B-A race, performed rather than argued.

        :class:`LoadedLibraries` fingerprints each source, then calls
        ``KineticsDatabase.load_libraries``, which re-opens that source **by
        pathname**. Those are two different reads of a mutable path, and the gap
        between them is a window. This test stands a writer in that window: it
        replaces the file, lets the loader consume the replacement, and puts the
        original back -- byte for byte, with the mtime forged back with ``utime`` --
        which is the shape that defeats a content hash entirely.

        The race is made deterministic instead of being run for real. A real
        concurrent writer would win this race only sometimes, and a test that
        depends on winning it is a flaky test that proves nothing on the runs where
        it loses; driving the swap from inside the loader's own call makes the
        window's *existence* the thing under test, which is what the finding is.

        Three assertions, in the order that makes the argument:

        1. the loader really did consume the replacement -- read off the object the
           loader produced, not inferred;
        2. **the sha256 of the path is identical before and after**, so the digest
           check that was the whole binding could not have seen it, and neither
           could the object-identity check: the object registered is the one this
           set loaded, because it is the one built from the impostor;
        3. the binding refuses anyway, on the inode identity, naming the interval.

        Everything is done on a COPY of the shipped library in ``tmp_path``. Nothing
        in ``RMG-database-plasma`` is written at any point.
        """
        scratch = tmp_path / 'libraries'
        shutil.copytree(os.path.join(self.libraries_path, RECOMBINATION),
                        str(scratch / RECOMBINATION))
        source = str(scratch / RECOMBINATION / 'reactions.py')
        pristine = str(tmp_path / 'pristine.py')
        shutil.copyfile(source, pristine)

        marker = 'SWAPPED-IN-THE-WINDOW'
        with open(pristine) as handle:
            text = handle.read()
        assert 'shortDesc' in text, 'the impostor cannot be built: no shortDesc to move'
        impostor_text = re.sub(r'^shortDesc = u".*?"$', 'shortDesc = u"{0}"'.format(marker),
                               text, count=1, flags=re.M | re.S)
        assert impostor_text != text, 'the impostor is byte-identical to the original'
        impostor = str(tmp_path / 'impostor.py')
        with open(impostor, 'w') as handle:
            handle.write(impostor_text)

        before = _file_fingerprint(source)
        stat_before = os.stat(source)
        real_load_libraries = KineticsDatabase.load_libraries
        consumed = {}

        def racing_load(self, path, libraries=None):
            """The writer, standing exactly in the window."""
            shutil.copyfile(impostor, source)
            try:
                result = real_load_libraries(self, path, libraries=libraries)
            finally:
                shutil.copyfile(pristine, source)
                os.utime(source, ns=(stat_before.st_atime_ns, stat_before.st_mtime_ns))
            consumed['library'] = self.libraries.get(RECOMBINATION)
            return result

        monkeypatch.setattr(KineticsDatabase, 'load_libraries', racing_load)
        with pytest.raises(AssertionError) as raised:
            LoadedLibraries(str(scratch), [RECOMBINATION])
        monkeypatch.undo()

        assert consumed.get('library') is not None, (
            'the loader never ran, so this test measured nothing')
        assert marker in (consumed['library'].short_desc or ''), (
            'the loader did not consume the replacement, so the race was not won and '
            'the refusal below is about something else. It read: {0!r}'.format(
                consumed['library'].short_desc))

        after = _file_fingerprint(source)
        assert after[0] == before[0], (
            'the premise moved: the file was NOT restored byte-for-byte, so a plain '
            'content hash would have caught this and there is no window to close')
        assert os.stat(source).st_mtime_ns == stat_before.st_mtime_ns, (
            'the premise moved: the mtime was not forged back, so an mtime check '
            'would have caught this')
        assert after[1] != before[1], (
            'the inode identity is unchanged across a replace-and-restore, so this '
            'repair has nothing to read and the window is not closable this way')

        message = str(raised.value)
        assert 'replaced and put back' in message, message
        assert 'loading them' in message, (
            'the refusal does not say WHEN the file moved; "while loading" and '
            '"since it was loaded" are different failures:\n{0}'.format(message))
        assert source in message, 'the refusal does not name the file'

        # The control, and it is load-bearing: the same construction with no writer
        # in the window must succeed. Without it this check would pass just as well
        # against a binding that refused every load.
        assert LoadedLibraries(str(scratch), [RECOMBINATION]).library(RECOMBINATION) is not None

    def test_every_entry_answers_to_the_index_it_carries(self):
        """The handle's own integrity, measured on every plasma library shipped.

        ``library.entries`` is keyed on the authored index and each entry also
        carries it on ``entry.index``. Selection reads the key; every assertion
        that follows reads the entry. If the two ever disagreed, the handle would
        name one entry and the assertions would describe another, which is the
        silent redirect this round exists to make impossible.

        Uniqueness is not asserted here because it cannot fail:
        ``KineticsLibrary.load_entry`` refuses a second entry under an index
        already taken (``library.py:680``). What that assert does NOT survive is
        ``python -O``, which strips it -- so this check names the reliance rather
        than leaving it implicit.
        """
        libraries = list(self.SHIPPED_CHANNELS) + list(self.LIBRARIES_NOT_SWEPT)
        database = LoadedLibraries(self.libraries_path, libraries)
        checked = 0
        for label in libraries:
            entries = database.libraries[label].entries
            assert entries, '{0} loaded no entries, so this check is vacuous on it'.format(label)
            for key, entry in entries.items():
                assert key == entry.index, (
                    'in library {0!r} the entry filed under key {1!r} carries index '
                    '{2!r} on itself. The selection handle and the entry disagree, so '
                    'a check naming {1!r} would assert against a different entry than '
                    'the one it selected.'.format(label, key, entry.index))
                assert isinstance(key, int), (
                    'in library {0!r} the entry key {1!r} is a {2}, not an int. The '
                    'index is an authored integer; a non-integer key means the loader '
                    'accepted something this handle cannot rely '
                    'on.'.format(label, key, type(key).__name__))
                checked += 1
        assert checked >= len(libraries), (
            'only {0} entr(ies) were checked across {1} librar(ies)'.format(
                checked, len(libraries)))

    #: The channels this file names, per library. The sweep's denominator is built
    #: from these **by provenance**, never from a frozen count.
    #:
    #: **A frozen count was tried first and the growth experiment refuted it.**
    #: The obvious repair to ``assert measured and set(values) == {2}`` -- which is
    #: satisfied by ONE entry as happily as by three, so a library that lost two of
    #: its three would pass the sweep it was supposed to fail -- is
    #: ``{IONISATION: 1, RECOMBINATION: 2}``. Measured: that form turns BOTH sweeps
    #: red against a database grown by one unrelated entry
    #: (``evidence/GROWTH-repaired.*``), which is the count-fragility this whole
    #: file exists to remove, reintroduced one level up. A census legitimately has
    #: a denominator; the denominator has to be *derived from the collection being
    #: censused*, and the claim has to be about coverage and named presence, not
    #: about size. Shrinkage is then caught by name (the missing channel is named,
    #: which a count could not do) and growth is tolerated.
    SHIPPED_CHANNELS = {IONISATION: (LITHIUM_IONISATION,),
                        RECOMBINATION: (LITHIUM_RECOMBINATION, ARGON_RECOMBINATION)}

    #: Why *this* sweep is over two libraries and not the five the database ships.
    #: ``PlasmaAir``, ``PlasmaAlkali`` and ``PlasmaArgon`` also carry plasma
    #: kinetics; they are outside this fixture because ``_database()`` loads
    #: exactly the two libraries this file's subject lives in, and widening the
    #: load would change what every other check in this class is measured against.
    #:
    #: **The reason given here until this round was wrong, and the measurement is
    #: what corrected it.** It said those three "are written in
    #: ``ElectronCollisionPlasma``". Measured on database ``96f2afa4a``: they are
    #: not. ``PlasmaAir`` holds 89 entries of which 3 are
    #: ``ElectronCollisionPlasma`` (the rest are ``Arrhenius``,
    #: ``TwoTemperaturePlasma`` and ``ThirdBody``); ``PlasmaAlkali`` holds 66 of
    #: which 5 are; only ``PlasmaArgon``, with its single entry, is written wholly
    #: in that class. Those libraries are also not all second order -- they carry
    #: unimolecular and three-body entries -- so a blanket "every entry is second
    #: order" sweep over them would be false, and that, not the class, is the real
    #: reason they cannot simply join this one.
    #:
    #: What they DO hold that this file needs exercising is the 9 real
    #: ``ElectronCollisionPlasma`` entries, and those are now swept, on their own
    #: terms, by :class:`TestTheShippedCrossSectionsAreRead`.
    LIBRARIES_NOT_SWEPT = ('PlasmaAir', 'PlasmaAlkali', 'PlasmaArgon')

    def _shipped_entries(self, database):
        """Every library reaction in the two loaded libraries, keyed by label, with
        the sweep's denominator established before anything is measured over it.

        Three statements, in this order, and each of them is a different failure:

        1. **Non-empty.** A library that loaded nothing makes every ``set(...)``
           assertion below vacuously true.
        2. **Every channel this file names is present**, found by provenance
           through :func:`select_library_reaction`, which raises naming every
           entry label the library does hold. This is what catches shrinkage, and
           it catches it *by name* -- "the argon recombination is gone" rather
           than "2 != 3".
        3. **Coverage**: every reaction the library yielded lands in the returned
           dict. Derived from the collection, not frozen, so growth is tolerated;
           a collision (two entries rendering identically) is caught, because that
           is how a sweep silently ends up over fewer entries than it loaded.
        """
        found = {}
        for label in (IONISATION, RECOMBINATION):
            reactions = database.libraries[label].get_library_reactions()
            assert reactions, (
                '{0} yielded no reactions at all, so every sweep over it is '
                'vacuously true'.format(label))
            for reference in self.SHIPPED_CHANNELS[label]:
                select_library_reaction(database.libraries[label], reference, reactions)
            before = len(found)
            for reaction in reactions:
                found['{0}: {1}'.format(label, reaction)] = reaction
            assert len(found) - before == len(reactions), (
                '{0} yielded {1} reactions but only {2} distinct keys, so two of '
                'them render identically and the sweep is silently over fewer '
                'entries than it loaded'.format(label, len(reactions),
                                                len(found) - before))
        return found

    def test_every_shipped_plasma_entry_is_second_order(self):
        """The descriptions' ``order=2`` is measured here, not assumed there.

        Three entries ship across the two libraries -- one ionisation, two
        recombinations -- and all three are second order: one heavy species and
        one electron, ``cm^3/(molecule*s)``. If a shipped rate is ever restated at
        another order this check names it, so the module-level descriptions cannot
        quietly go stale without anything saying so.

        The denominator is established by :meth:`_shipped_entries`, not observed.
        The previous form of this check was ``assert measured and set(...) == {2}``,
        which is green on any non-empty subset -- so the library could have lost
        two of its three entries and the sweep would still have reported that
        every entry it found was second order. That sentence is true and useless.
        """
        database = self._database()
        measured = {key: rate_order(reaction.kinetics)
                    for key, reaction in self._shipped_entries(database).items()}
        assert set(measured.values()) == {2}, measured

    def test_every_shipped_plasma_entry_has_the_response_its_channel_implies(self):
        """The same sweep on the new key, with the same asserted denominator.

        Each entry's Te response must be the one its chemistry implies: the
        ionisation entry rises (a threshold process), both recombination entries
        fall (a faster electron is less likely to be captured). This is what makes
        the ``response=`` on the module-level descriptions a measurement of the
        shipped fits rather than a label copied off their names.
        """
        database = self._database()
        expected = {IONISATION: RISES, RECOMBINATION: FALLS}
        measured = {key: electron_temperature_response(reaction.kinetics)
                    for key, reaction in self._shipped_entries(database).items()}
        for key, response in measured.items():
            label = key.split(':')[0]
            assert response == expected[label], (
                '{0} does not have the Te response its channel implies: wanted {1}, '
                'measured {2}\n  all: {3}'.format(key, expected[label], response,
                                                  measured))

    def test_the_lithium_recombination_is_radiative_not_an_ionisation_rate(self):
        """The mutation this round exists for, measured against the real fits.

        Both objects below are shipped by this database for this element:
        ``BadnellRRArrhenius(Z=3, N=2)`` is the radiative recombination rate the
        entry carries, and ``VoronovEIArrhenius(Z=3, N=3)`` is the ionisation rate
        the other library carries. Pasting the second onto the recombination entry
        -- with ``electrons = -1`` so even the bookkeeping is untouched -- is the
        substitution closest to the accepted case that this database makes
        possible: nothing is invented, one shipped fit stands in for another.

        The first block is the premise: every key except the response accepts the
        substitution. The second is the discriminator. The third puts the
        substitute *on the shipped entry* and asserts that the entry is still
        selected -- it is still that entry, and pretending otherwise was the old
        design's mistake -- and that the assertion refuses it by name.
        """
        database = self._database()
        radiative = self._select(database, LITHIUM_RECOMBINATION)
        substitute = VoronovEIArrhenius(Z=3, N=3, electrons=-1)

        # The premise: every key but one accepts the substitution.
        assert rate_order(substitute) == rate_order(radiative.kinetics) == 2, (
            'the rate order already separates these, so the response adds nothing')
        # ``KineticsModel.electrons`` is a ``ScalarQuantity``, not an int, and
        # ``ScalarQuantity`` has no value equality: measured, ``q == -1`` is False
        # and so is ``q1 == q2`` for two quantities of the same value. Comparing
        # the objects would make this premise line vacuously "different" and the
        # premise would read as proven when it had not been asked. ``.value`` is
        # the number.
        assert substitute.electrons.value == radiative.electrons == -1, (
            'the net electron count already separates these')
        assert (getattr(substitute.A, 'units', None)
                == getattr(radiative.kinetics.A, 'units', None)), (
            'the A-factor units already separate these')

        # The discriminator.
        assert electron_temperature_response(radiative.kinetics) == FALLS
        assert electron_temperature_response(substitute) == RISES

        # The shipped entry passes; the substitute, on the same entry, is refused
        # by name rather than sending the selection elsewhere.
        LITHIUM_RECOMBINATION.assert_as_shipped(radiative)
        with pytest.raises(AssertionError) as raised:
            LITHIUM_RECOMBINATION.assert_as_shipped(radiative, kinetics=substitute)
        message = str(raised.value)
        assert 'wrong KINETICS CLASS' in message or 'wrong Te RESPONSE' in message, message

    def test_the_lithium_recombination_is_the_radiative_channel_not_the_three_body_one(self):
        """The entry names a *channel*, and the database is asked which one.

        Both channels are written ``[Lip] => [Li]`` with ``electrons = -1``, so no
        participant-level key can tell them apart. The shipped entry is the
        second-order radiative one; the third-order collisional one is not in this
        database at all (``PlasmaRadiativeRecombination``'s ``longDesc`` says so:
        "no shipped fit anywhere in this database").

        Note what this check can no longer say, and why that is an improvement.
        It used to assert that *asking for* a third-order lithium recombination
        found nothing -- but "nothing matched" is also what a renamed entry, a
        failed load or an empty library produce, so the old form was green for
        several reasons that had nothing to do with the channel. The entry is now
        named directly and its order asserted, which is a statement about the
        shipped data rather than about a search.
        """
        database = self._database()
        radiative = self._select(database, LITHIUM_RECOMBINATION)
        assert rate_order(radiative.kinetics) == 2, (
            'the shipped lithium recombination is no longer second order, so it is no '
            'longer the radiative channel: {0!r}'.format(radiative.kinetics))
        assert electron_temperature_response(radiative.kinetics) == FALLS
        assert radiative.kinetics.is_identical_to(BadnellRRArrhenius(Z=3, N=2)), (
            'the shipped entry is no longer the Badnell radiative fit for lithium')

    def test_the_two_channels_are_not_a_reverse_pair(self):
        """The measurement the whole repair turns on. Their net counts are equal
        and opposite -- the relation a reverse pair has -- but their placements
        are not mirror images, because the reverse of the ionisation is the
        three-body channel and not the radiative one."""
        database = self._database()
        ionisation = self._select(database, LITHIUM_IONISATION)
        recombination = self._select(database, LITHIUM_RECOMBINATION)

        assert ionisation.electrons == -recombination.electrons

        ion_counts = get_electron_placement_counts(ionisation)
        rec_counts = get_electron_placement_counts(recombination)
        assert ion_counts == (1, 2)
        assert rec_counts == (1, 0)
        assert (ion_counts[1], ion_counts[0]) != rec_counts

    @pytest.mark.parametrize('source_first', [True, False])
    def test_both_channels_enter_the_model(self, source_first):
        """Both orders, because the drop was order-dependent: whichever was
        offered second was the one that disappeared."""
        database = self._database()
        order = [IONISATION, RECOMBINATION] if source_first else [RECOMBINATION, IONISATION]
        verdicts = {}
        with _GlobalKineticsDatabase(database.kinetics_database):
            model = CoreEdgeReactionModel()
            model.kinetics_database = database.kinetics_database
            for label in order:
                reaction = self._select(database, LITHIUM_CHANNEL[label])
                _, is_new = model.make_new_reaction(reaction,
                                                    generate_thermo=False,
                                                    generate_kinetics=False)
                verdicts[label] = is_new

        assert verdicts[IONISATION] is True, 'the ionisation source did not enter the model'
        assert verdicts[RECOMBINATION] is True, 'the recombination sink did not enter the model'

    def test_the_cation_has_a_loss_channel_in_the_model(self):
        """The deliverable, stated as chemistry rather than as a count: the cation
        appears on the reactant side of some reaction in the model."""
        database = self._database()
        with _GlobalKineticsDatabase(database.kinetics_database):
            model = CoreEdgeReactionModel()
            model.kinetics_database = database.kinetics_database
            for label in (IONISATION, RECOMBINATION):
                model.make_new_reaction(self._select(database, LITHIUM_CHANNEL[label]),
                                        generate_thermo=False, generate_kinetics=False)

        produced = [product for rxn in model.new_reaction_list for product in rxn.products
                    if product.molecule and product.molecule[0].get_net_charge() > 0]
        assert produced, 'no cation is produced at all'

        unconsumed = cations_produced_but_never_consumed(model.new_reaction_list)
        assert not unconsumed, (
            'these cations are produced but never consumed: {0}'.format(
                sorted(_identify(species) for species in unconsumed)))

    def test_a_label_collision_does_not_satisfy_the_cation_loss_check(self):
        """The cation-loss check reads species, not the strings that name them.

        Until this round it collapsed both sides to ``species.label`` and compared
        those, which makes an unrelated reactant that happens to render under a
        produced cation's label sufficient. That is not a contrived collision in
        this format: a label here is a rendering of the chemistry, so a species
        rebuilt on the other side of any boundary carries the same one, and the
        check would have reported a loss channel that the model does not have.

        The fixture is the nearest wrong thing to the accepted case, not an
        obviously broken one. The model produces a cation and consumes NOTHING;
        a second reaction consumes a **different** species object carrying the
        **same label**. Under the label comparison that model passes. Under
        identity it does not, and the produced cation is named.

        The positive control is the same fixture with the real object consumed --
        without it this check would pass just as well against a predicate that
        refused everything.
        """
        cation = _structure(LITHIUM_CATION, label='Lip')
        neutral = _structure(LITHIUM_ATOM, label='Li')
        impostor = _structure(LITHIUM_ATOM, label='Lip')
        assert impostor is not cation and impostor.label == cation.label, (
            'the fixture is not the case this is about: the impostor must be a '
            'DIFFERENT object carrying the SAME label')

        ionisation = _reaction([neutral], [cation], electrons=1)
        red_herring = _reaction([impostor], [neutral], electrons=0)
        unconsumed = cations_produced_but_never_consumed([ionisation, red_herring])
        assert [species for species in unconsumed] == [cation], (
            'a different species sharing the cation\'s label counts as consuming it, '
            'so the check is about the string and not about the chemistry')

        genuine = _reaction([cation], [neutral], electrons=-1)
        assert cations_produced_but_never_consumed([ionisation, genuine]) == [], (
            'the control moved: consuming the cation itself no longer satisfies the '
            'check, so it now refuses everything')


#: The three anchored entries, paired with a rate law carrying exactly what each
#: one is asserted to carry. Built here rather than read from the database so the
#: adequacy checks below run without one: they are about the SHAPE of the anchor
#: set against a rate law's parameters, which is a property of this file's design
#: and not of any particular checkout.
#:
#: ``test_the_reference_rate_laws_are_what_the_database_ships`` ties them back to
#: the real entries, so this is a fixture and not a fourth transcription.
def _selected_entry_carrying(reference, kinetics):
    """``(reaction, entry)`` for a one-entry library holding `reference`'s entry.

    Built from the reference's own asserted attributes, so the fixture carries
    exactly what the entry is supposed to carry apart from the rate law under
    test. That is what makes a refusal attributable: anything
    :meth:`ShippedEntry.assert_as_shipped` reports on this fixture is about the
    kinetics, because nothing else was moved.

    Selection is by provenance -- library label and entry index -- and runs
    through the same ``get_library_reactions`` path the database-backed checks
    use, so this is not a second selector.
    """
    reaction = _reaction(
        [_structure(adjacency_list) for adjacency_list in reference.reactants],
        [_structure(adjacency_list) for adjacency_list in reference.products],
        electrons=reference.electrons, family=reference.library,
        reversible=reference.reversible, duplicate=reference.duplicate,
        allow_pdep_route=reference.allow_pdep_route,
        elementary_high_p=reference.elementary_high_p)
    library = _fixture_library(reference.library, [
        _fixture_entry(reference.entry_index, reference.entry_label,
                       reaction, kinetics)])
    return (select_library_reaction(library, reference),
            select_entry(library, reference))


def _argon_reference_law():
    return TwoTemperaturePlasma(
        A=(_SVS_ARGON['A_rad'], 'cm^3/(molecule*s)'), n=-_SVS_ARGON['eta'],
        Ea_g=(0.0, 'kJ/mol'), Ea_e=(0.0, 'kJ/mol'), T0=(1.0e4, 'K'),
        electrons=-1, Tmin=(1.0e4, 'K'), Tmax=(1.0e8, 'K'))


ANCHORED_ENTRIES = (
    (LITHIUM_IONISATION, lambda: VoronovEIArrhenius(Z=3, N=3)),
    (LITHIUM_RECOMBINATION, lambda: BadnellRRArrhenius(Z=3, N=2)),
    (ARGON_RECOMBINATION, _argon_reference_law),
)


def _anchored_entry_ids(case):
    return case[0].name.replace(' ', '-')


class TestTheAnchorSetDeterminesEveryParameter:
    """Is the anchor set ADEQUATE -- computed here, not claimed in a docstring.

    **Why this class exists at all.** The same defect has now occurred twice in
    this file, in the same shape: an evaluation set that lost a dimension, which
    nobody re-derived. Round 91's check collapsed the rate law onto ``T == Te``
    and lost the gas-temperature axis; its repair evaluated at two ``Te`` values
    at one ``Tgas`` and lost it again, so ``A`` and ``Ea_g`` traded against each
    other exactly and a mutation 44x wrong one axis over was identical to nine
    decimals where the anchors looked.

    Both times the anchor set was described as sufficient and was not. So
    sufficiency is no longer described: it is COMPUTED from the anchor grid and
    the entry's own rate law, and asserted. Removing an anchor point, flattening
    an axis, or adding a parameter to one of these classes is red by
    construction rather than red the next time somebody re-derives the algebra.

    **What "adequate" is taken to mean here.** Not "rank 4". Rank is necessary
    and not sufficient -- a formally full-rank design with a bad condition number
    still hides a perturbation below the tolerance in force -- and for two of the
    three rate laws there is no exact design matrix at all, because Badnell's and
    Voronov's forms are not linear in their parameters. So the property asserted
    is the one that survives both objections and is measurable for all three:

        every free parameter of the entry's own rate law is pinned by at least
        one mechanism this file actually runs, and for each parameter the
        smallest perturbation that mechanism refuses is measured and recorded.

    Two mechanisms do the pinning, and they are not redundant:

    * the ANCHORS, which catch a parameter through the rate it produces when put
      through production's own evaluator -- and so also catch a rate law whose
      fields are right and whose function of them is wrong;
    * the DIRECT PARAMETER ASSERTION (``ShippedEntry.rate_law_parameters``),
      which catches a parameter by reading it -- and so also catches the
      parameters the anchors provably cannot reach.

    ``test_every_rate_law_parameter_is_pinned_by_some_mechanism`` is the
    assertion that the union covers everything.
    """

    @pytest.mark.parametrize('case', ANCHORED_ENTRIES, ids=_anchored_entry_ids)
    def test_the_parameter_census_matches_each_rate_laws_own_arity(self, case):
        """:data:`RATE_LAW_PARAMETERS` is a hand-written census, so it is tied back.

        A table in this file enumerating another module's fields is the mirror
        shape this campaign keeps finding stale: N normative copies compared to
        each other and none to the implementation. The tie-back available here is
        the class's own ``__reduce__`` tuple, which pickling forces to enumerate
        every field that matters.

        **The arity alone was not enough, and this round says why.** The check used
        to pin the total arity plus "every recorded parameter slot is in range and
        readable". A bookkeeping slot replaced by a new functional parameter passes
        both: the count does not move and no recorded slot becomes invalid. So the
        census now accounts for EVERY slot and asserts what each one holds --
        see :func:`assert_reduce_slots_match`.
        """
        reference, build = case
        assert_reduce_slots_match(build())

    def test_a_bookkeeping_slot_that_became_a_parameter_is_refused(self):
        """The mutation the arity check could not see, run through the real check.

        The three shipped rate laws are compiled ``cdef`` classes, so a genuinely
        changed slot cannot be produced without editing and rebuilding
        ``rmgpy/kinetics/arrhenius.pyx`` -- which this branch does not touch. The
        subject is therefore a stand-in whose ``__reduce__`` has the same shape as
        Voronov's: five parameters, ``electrons``, three slots the class hard-codes
        as ``None``, and the bookkeeping tail.

        The mutation is the nearest wrong thing to the accepted case, not an
        obviously broken one: **one** of those hard-coded ``None`` slots now carries
        a value, and nothing else moves. The arity is identical, every recorded
        parameter slot is still in range, and every parameter still reads back --
        so the check as it stood was green on it. It is put through
        :func:`assert_reduce_slots_match`, the same function the three real rate
        laws go through, and not through a second copy of the rule.

        The control is the unmutated stand-in: without it this would pass equally
        against a check that refused everything.
        """

        class _StandIn:
            """A rate law with Voronov's reduce shape."""

            def __init__(self, reserved=None):
                self.A, self.P, self.X, self.K, self.dE = 1.0, 0.0, 0.4, 0.4, 5.4
                self.electrons = 1
                self.Tmin, self.Tmax, self.Pmin, self.Pmax = 10.0, 1.0e8, None, None
                self.uncertainty, self.solute, self.comment = None, None, ''
                self._reserved = reserved

            def __reduce__(self):
                return (_StandIn, (self.A, self.P, self.X, self.K, self.dE,
                                   self.electrons,
                                   self._reserved, None, None,
                                   self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                                   self.uncertainty, self.solute, self.comment))

        census = dict(RATE_LAW_PARAMETERS['VoronovEIArrhenius'])
        assert_reduce_slots_match(_StandIn(), census=census)

        mutated = _StandIn(reserved=0.37)
        assert len(mutated.__reduce__()[1]) == len(_StandIn().__reduce__()[1]), (
            'the mutation changed the arity, so the old check would have caught it '
            'and this is not the case under test')
        for _, slot in census['parameters']:
            assert slot < len(census['reduce_slots']), (
                'the mutation invalidated a recorded parameter slot, so the old '
                'check would have caught it')
            assert getattr(mutated, census['reduce_slots'][slot], None) is not None, (
                'the mutation made a recorded parameter unreadable, so the old check '
                'would have caught it')

        with pytest.raises(AssertionError) as raised:
            assert_reduce_slots_match(mutated, census=census)
        message = str(raised.value)
        assert 'slot 6' in message, (
            'the refusal does not name which slot moved:\n{0}'.format(message))
        assert '0.37' in message, (
            'the refusal does not show what appeared there:\n{0}'.format(message))

    @pytest.mark.parametrize('case', ANCHORED_ENTRIES, ids=_anchored_entry_ids)
    def test_the_slot_map_and_the_parameter_list_agree(self, case):
        """The two halves of the census are checked against each other as well.

        ``parameters`` is what ``perturb_rate_law_parameter`` rebuilds the class
        from, and ``reduce_slots`` is what the tie-back reads. A disagreement
        between them would mean every measured sensitivity for the mis-pointed
        parameter was a measurement of a different field, silently.
        """
        reference, build = case
        kinetics = build()
        census = RATE_LAW_PARAMETERS[type(kinetics).__name__]
        assert tuple(name for name, _ in census['parameters']) \
            == rate_law_parameter_names(kinetics)
        for name, slot in census['parameters']:
            assert census['reduce_slots'][slot] == name

    @pytest.mark.parametrize('case', ANCHORED_ENTRIES, ids=_anchored_entry_ids)
    def test_the_anchor_grid_has_the_axes_the_rate_law_needs(self, case):
        """Three distinct Te always; two distinct Tgas exactly when they buy something.

        This is the structural half of the repair -- the part that makes losing a
        dimension red rather than invisible. It asserts the axes of the grid the
        entry is actually anchored on, and that the gas-temperature axis is
        present on precisely the entries whose rate law varies with the gas
        temperature.

        The 2x2 grid is the trap this guards against and it deserves naming: the
        Kossyi form is additive in ``ln Te`` and ``1/T`` with no interaction term,
        so a tensor product of two values on each axis is rank-deficient however
        far apart the two values are. Two axes are not enough; three points on one
        of them is what makes the difference.
        """
        reference, build = case
        kinetics = build()
        electron_temperatures = {anchor.temperature
                                 for anchor in reference.rate_anchors}
        gas_temperatures = {anchor.gas_temperature
                            for anchor in reference.rate_anchors}
        wanted = len(rate_law_parameter_names(kinetics))
        assert len(electron_temperatures) >= min(wanted, 6), (
            '{0} is anchored at {1} distinct electron temperature(s) and its rate '
            'law has {2} free parameters. An anchor set can separate at most as '
            'many independent directions as it has points, so fewer points than '
            'parameters is under-determined before conditioning is even '
            'considered.'.format(reference.name, len(electron_temperatures), wanted))
        needs_gas_axis = rate_law_admits_gas_temperature_dependence(kinetics)
        assert (len(gas_temperatures) >= 2) == needs_gas_axis, (
            '{0}\'s rate law {1} a function of the gas temperature (measured), and '
            'it is anchored at {2} distinct gas temperature(s). A rate law that '
            'varies with Tgas and is anchored at one of them cannot have A and '
            'Ea_g separated at all; a rate law that does not vary with Tgas gains '
            'nothing from a second and should not claim to.'.format(
                reference.name, 'IS' if needs_gas_axis else 'is NOT',
                len(gas_temperatures)))
        assert set(reference.rate_anchors and
                   [(a.gas_temperature, a.temperature)
                    for a in reference.rate_anchors]) == set(anchor_grid(kinetics)), (
            '{0}\'s anchors are not the grid anchor_grid() computes for its own '
            'rate law, so the explicit list in the entry definition has drifted '
            'from the design that justifies it.'.format(reference.name))

    @pytest.mark.parametrize('case', ANCHORED_ENTRIES, ids=_anchored_entry_ids)
    def test_every_rate_law_parameter_is_pinned_by_some_mechanism(self, case):
        """The sufficiency property: every parameter is READ, not merely evaluated.

        This is the half the anchors cannot supply however good the grid is, and
        the reason is structural rather than a matter of tolerance. The anchors
        constrain the rate law only through the values it produces, so any
        perturbation lying in the null space of the sensitivity matrix is
        invisible to them at every point, at every tolerance, for ever. For the
        argon entry that null space is not empty and never can be: ``A`` and
        ``T0`` enter the Kossyi form only through ``ln A - n*ln T0``, so they are
        a reparametrisation of one another.

        A direct assertion has no null space, because it reads each field. So the
        requirement asserted here is total coverage: every free parameter of every
        anchored entry's rate law has a value recorded for it in
        ``rate_law_parameters``. Not "most", and not "the ones the anchors miss" --
        which would need the null space to be re-derived by hand every time the
        grid moved, and re-deriving by hand is what failed twice already.
        """
        reference, build = case
        kinetics = build()
        names = set(rate_law_parameter_names(kinetics))
        asserted = set(reference.rate_law_parameters)
        assert names <= asserted, (
            '{0}: {1} has free parameter(s) {2} that no direct assertion reads. '
            'The anchors alone cannot close this: any perturbation in the null '
            'space of their sensitivity matrix is invisible to them at every '
            'point and every tolerance, and for this rate law that null space is '
            'not empty.'.format(reference.name, type(kinetics).__name__,
                                ', '.join(sorted(names - asserted))))
        assert asserted <= names, (
            '{0}: a value is asserted for {1}, which is not a free parameter of '
            '{2} -- so nothing reads it and the assertion is dead.'.format(
                reference.name, ', '.join(sorted(asserted - names)),
                type(kinetics).__name__))

    @pytest.mark.parametrize('case', ANCHORED_ENTRIES, ids=_anchored_entry_ids)
    def test_the_anchors_separate_every_parameter_they_are_not_declared_blind_to(self, case):
        """The anchors' own adequacy: full column rank, minus the declared degeneracies.

        The anchors' job is not to pin the numbers -- the direct assertion above
        does that, and does it without a null space. Their job is to check the
        FORM and the DISPATCH: that these fields, put through production's own
        evaluator, reproduce the published rate. A field assertion is blind to a
        rate law that carries the right parameters and computes the wrong function
        of them, or that production routes to the wrong evaluator, and those are
        real failure modes this campaign has met.

        For the anchors to do that job they have to remain sensitive to each
        parameter separately, and that is what is asserted here: the rank of the
        sensitivity matrix over the anchor grid equals the number of free
        parameters, less exactly the degeneracies the entry DECLARES. Flattening
        an axis or dropping a point costs rank and turns this red -- which is the
        structural property the last two rounds of this defect lacked.

        A declared degeneracy has to be an exact algebraic identity, not an
        observation that two columns look similar, and the null direction is
        checked to lie inside the declared set rather than merely to exist.
        """
        reference, build = case
        kinetics = build()
        names, grid, jacobian = anchor_sensitivity(kinetics)
        degeneracies = reference.anchor_degeneracies
        expected_rank = len(names) - sum(len(group) - 1 for group in degeneracies)
        singular = np.linalg.svd(jacobian, compute_uv=False)
        rank = int(np.sum(singular > 1.0e-8 * singular[0]))
        assert rank == expected_rank, (
            '{0}: the anchors resolve {1} independent directions of {2}\'s {3} '
            'parameters, and {4} was expected ({3} parameters less the declared '
            'degeneracies {5}). Singular values {6}. If an axis was flattened or a '
            'point dropped, restore it; if a new exact degeneracy has been '
            'introduced, it has to be declared and covered.'.format(
                reference.name, rank, type(kinetics).__name__, len(names),
                expected_rank, degeneracies or '(none)',
                np.array2string(singular, precision=3)))
        if degeneracies:
            declared = {name for group in degeneracies for name in group}
            _, _, right = np.linalg.svd(jacobian)
            for direction in right[rank:]:
                carried = {name for name, weight in zip(names, direction)
                           if abs(weight) > 1.0e-6}
                assert carried <= declared, (
                    '{0}: the anchors cannot separate {1}, which is outside the '
                    'declared degeneracies {2}. A parameter combination that is '
                    'invisible to every anchor is a quantity this file cannot '
                    'see; work out whether it is an exact identity (declare it) '
                    'or a conditioning failure (fix the grid).'.format(
                        reference.name, ', '.join(sorted(carried)), degeneracies))

    @pytest.mark.parametrize('case', ANCHORED_ENTRIES, ids=_anchored_entry_ids)
    def test_the_measured_per_parameter_resolution_is_recorded(self, case):
        """What the anchors actually resolve, printed and bounded.

        The brief for this round asked for the smallest relative perturbation the
        design catches PER PARAMETER at the tolerance in force, rather than a
        claim that the design is adequate. That number is
        ``ANCHOR_RTOL / max_i abs(d ln k / d delta)`` over the anchor grid, and it
        is asserted against a bound so that a future edit which formally keeps the
        rank but wrecks the conditioning is red.

        Measured on the grid this file now asserts, with ``ANCHOR_RTOL = 1e-5``
        (``evidence/anchor_adequacy.py`` re-runs all of it):

            argon    A 1.0e-5  n 2.5e-6  Ea_g 2.5e-5 (kJ/mol)  Ea_e 2.5e-5
                     T0 1.5e-5 -- but see the A/T0 degeneracy: a LONE T0 move is
                     seen at 1.5e-5 and the coordinated A/T0 trade never is
            Badnell  A 1.0e-5  B 6.0e-6  T0 1.2e-5  T1 2.9e-5  C 1.7e-5  T2 4.7e-5
            Voronov  A 1.0e-5  P 5.6e-6  X 1.0e-5   K 5.6e-6   dE 2.8e-6

        The bound below is 1e-3, more than an order looser than the worst of
        those, so it is a tripwire against a design that stops resolving -- not a
        restatement of today's numbers, which would go red on float noise.
        """
        reference, build = case
        kinetics = build()
        thresholds = smallest_detected_perturbations(kinetics)
        reachable = {name: value for name, value in thresholds.items()
                     if value is not None}
        assert reachable, (
            '{0}: the anchors resolve none of its parameters at all'.format(
                reference.name))
        worst_name = max(reachable, key=reachable.get)
        assert reachable[worst_name] <= 1.0e-3, (
            '{0}: the anchors now need a perturbation of {1:.3g} in {2} before they '
            'refuse it, which is {3:.0f}x looser than the worst measured when this '
            'bound was set. The grid is formally still a grid and has stopped '
            'resolving; re-derive it rather than loosening this '
            'bound.'.format(reference.name, reachable[worst_name], worst_name,
                            reachable[worst_name] / 1.0e-3))

    @pytest.mark.parametrize('case', ANCHORED_ENTRIES, ids=_anchored_entry_ids)
    def test_each_parameter_perturbed_alone_is_refused(self, case):
        """Every parameter, moved BY ITSELF, and shown caught end to end.

        A combined mutation proves less than N separate ones, because any single
        parameter in it could be carrying the whole signal while the others are
        free. So each parameter is moved on its own, through the real
        :meth:`ShippedEntry.assert_as_shipped` on a real selected entry, and the
        refusal is required to name the parameter or the anchor that saw it.

        The perturbation is 1% -- far above every threshold recorded in
        ``test_the_measured_per_parameter_resolution_is_recorded`` and far above
        the 1e-9 the direct assertion uses, so a failure here is a hole and not a
        tolerance argument.
        """
        reference, build = case
        kinetics = build()
        control_reaction, control_entry = _selected_entry_carrying(reference, kinetics)
        reference.assert_as_shipped(control_reaction, entry=control_entry)

        for name in rate_law_parameter_names(kinetics):
            moved = perturb_rate_law_parameter(kinetics, name, 0.01)
            reaction, entry = _selected_entry_carrying(reference, moved)
            with pytest.raises(AssertionError) as raised:
                reference.assert_as_shipped(reaction, entry=entry)
            message = str(raised.value)
            assert ('wrong RATE MAGNITUDE' in message
                    or 'wrong RATE LAW PARAMETER' in message
                    or 'wrong FIT' in message), (
                '{0}: moving {1} by 1% is not refused as a magnitude, a parameter '
                'or a fit -- it was caught for some other reason, or not at '
                'all:\n{2}'.format(reference.name, name, message))
            assert name in message, (
                '{0}: moving {1} by 1% is refused, but the message never names '
                '{1}, so a reader of a red run cannot tell which parameter '
                'moved:\n{2}'.format(reference.name, name, message))

    def test_the_gas_temperature_trade_between_a_and_ea_g_is_refused(self):
        """The review's mutation: ``Ea_g = 10 kJ/mol`` with ``A`` rescaled to match.

        This is the one the old anchor set could not see, and the reason this
        round happened. Both parameters move together, chosen so the product
        ``A * exp(-Ea_g/(R*SOLVER_TGAS))`` is unchanged -- which makes the
        mutation exactly identical to the shipped rate everywhere on the line
        ``Tgas = SOLVER_TGAS``, where every anchor used to live, and 17x to 44x
        wrong off it.

        The check asserts BOTH halves: that the mutation really is invisible at
        the old anchors (otherwise it is not this mutation), and that the new grid
        refuses it.
        """
        shipped = _argon_reference_law()
        rescaled = 10.0e3    # J/mol
        traded = TwoTemperaturePlasma(
            A=(_SVS_ARGON['A_rad'] / math.exp(-rescaled / (constants.R * SOLVER_TGAS)),
               'cm^3/(molecule*s)'),
            n=-_SVS_ARGON['eta'], Ea_g=(10.0, 'kJ/mol'), Ea_e=(0.0, 'kJ/mol'),
            T0=(1.0e4, 'K'), electrons=-1, Tmin=(1.0e4, 'K'), Tmax=(1.0e8, 'K'))

        # Half one: it really is the mutation the old design could not see.
        for temperature in ANCHOR_TEMPERATURES:
            good = evaluate_as_the_solver_does(shipped, temperature, SOLVER_TGAS)
            bad = evaluate_as_the_solver_does(traded, temperature, SOLVER_TGAS)
            assert abs(bad - good) <= 1e-9 * good, (
                'the trade is visible at Tgas = {0} K, Te = {1} K, so this fixture '
                'is no longer the mutation that motivated the repair'.format(
                    SOLVER_TGAS, temperature))
        assert electron_temperature_response(traded) == FALLS
        assert rate_order(traded) == rate_order(shipped) == 2

        # Half two: off that line it is grossly wrong, and the new grid is there.
        off_axis = ANCHOR_GAS_TEMPERATURES[1]
        ratio = (evaluate_as_the_solver_does(traded, ANCHOR_TEMPERATURES[0], off_axis)
                 / evaluate_as_the_solver_does(shipped, ANCHOR_TEMPERATURES[0], off_axis))
        assert ratio > 10.0, (
            'the trade is only {0:.3g}x wrong at the second gas temperature, so '
            'that axis is not buying what this test claims'.format(ratio))

        reaction, entry = _selected_entry_carrying(ARGON_RECOMBINATION, traded)
        with pytest.raises(AssertionError) as raised:
            ARGON_RECOMBINATION.assert_as_shipped(reaction, entry=entry)
        message = str(raised.value)
        assert 'Tgas = {0} K'.format(off_axis) in message, (
            'the trade is refused, but not by an anchor at the second gas '
            'temperature -- so the axis this round added is not what caught '
            'it:\n{0}'.format(message))
        assert 'wrong RATE LAW PARAMETER' in message and 'Ea_g' in message, (
            'the trade is refused by the anchors but the parameter assertion '
            'never names Ea_g:\n{0}'.format(message))

    def test_the_combined_slope_and_electron_activation_mutation_is_refused(self):
        """The review's second mutation: ``n``, ``Ea_e`` and ``A`` moved together.

        A mutation chosen to be mild everywhere the old checks looked and 20.5%
        wrong at Te = 1e4 K. It is kept as a named case because a combined
        mutation is what a plausible bad edit looks like -- three fields adjusted
        to keep something else fixed -- and because it exercises the anchors and
        the parameter assertion at once.
        """
        shipped = _argon_reference_law()
        combined = TwoTemperaturePlasma(
            A=(_SVS_ARGON['A_rad'] * 1.05, 'cm^3/(molecule*s)'),
            n=-_SVS_ARGON['eta'] + 0.02, Ea_g=(0.0, 'kJ/mol'),
            Ea_e=(0.5, 'kJ/mol'), T0=(1.0e4, 'K'), electrons=-1,
            Tmin=(1.0e4, 'K'), Tmax=(1.0e8, 'K'))
        probe = 1.0e4
        deviation = abs(
            evaluate_as_the_solver_does(combined, probe, SOLVER_TGAS)
            / evaluate_as_the_solver_does(shipped, probe, SOLVER_TGAS) - 1.0)
        assert deviation > 0.10, (
            'the combined mutation is only {0:.1%} wrong at Te = {1} K, so it is no '
            'longer the case the review named'.format(deviation, probe))

        reaction, entry = _selected_entry_carrying(ARGON_RECOMBINATION, combined)
        with pytest.raises(AssertionError) as raised:
            ARGON_RECOMBINATION.assert_as_shipped(reaction, entry=entry)
        message = str(raised.value)
        assert 'wrong RATE MAGNITUDE' in message, (
            'the combined mutation is not refused on magnitude:\n{0}'.format(message))
        for name in ('A', 'n', 'Ea_e'):
            assert name in message, (
                'the refusal never names {0}, one of the three fields that '
                'moved:\n{1}'.format(name, message))

    def test_a_parameter_no_anchor_can_reach_is_pinned_by_a_direct_assertion(self):
        """``T0``: the parameter the brief for this round did not know was free.

        The repair was scoped as four free parameters -- ``A``, ``n``, ``Ea_g``,
        ``Ea_e``. There are five. ``T0`` enters the Kossyi form only through
        ``ln A - n*ln T0``, so it and ``A`` are a reparametrisation of one
        another: no anchor set, at any tolerance, on any number of points, at any
        pair of temperatures, can separate them. Measured, and asserted here:
        moving ``T0`` from 1e4 K to 1 K with ``A`` multiplied by ``(1/1e4)^n``
        leaves every anchor identical.

        This is why sufficiency had to be computed rather than argued. A
        rank-based argument over a design matrix would have reported rank 4 and
        called the job done; the parameter it cannot see does not appear in the
        rank at all unless the parameter census is complete and checked against
        it. The anchors are not the pin for ``T0`` and cannot be made into one;
        the direct parameter assertion is.
        """
        shipped = _argon_reference_law()
        moved_t0 = 1.0
        traded = TwoTemperaturePlasma(
            A=(_SVS_ARGON['A_rad'] * (moved_t0 / 1.0e4) ** -_SVS_ARGON['eta'],
               'cm^3/(molecule*s)'),
            n=-_SVS_ARGON['eta'], Ea_g=(0.0, 'kJ/mol'), Ea_e=(0.0, 'kJ/mol'),
            T0=(moved_t0, 'K'), electrons=-1, Tmin=(1.0e4, 'K'), Tmax=(1.0e8, 'K'))

        # The anchors cannot see it, and that is a property of the algebra.
        for gas, electron in anchor_grid(shipped):
            good = evaluate_as_the_solver_does(shipped, electron, gas)
            bad = evaluate_as_the_solver_does(traded, electron, gas)
            assert abs(bad - good) <= 1e-11 * good, (
                'the T0/A trade is visible at (Tgas={0}, Te={1}), so it is not the '
                'exact reparametrisation this test is about'.format(gas, electron))
        # And it is the PAIR that is invisible, not T0 by itself. A lone T0 move
        # is seen (measured: 1.5e-5 is the smallest the anchors refuse), which is
        # exactly why a per-parameter reading of the sensitivity matrix reports
        # T0 as covered and is the wrong instrument here. The null space is.
        thresholds = smallest_detected_perturbations(shipped)
        assert thresholds['T0'] is not None, (
            'a lone T0 perturbation has become invisible to the anchors, which is '
            'a different and larger defect than the trade this test is about')
        names, _, jacobian = anchor_sensitivity(shipped)
        _, singular, right = np.linalg.svd(jacobian)
        rank = int(np.sum(singular > 1.0e-8 * singular[0]))
        assert rank == len(names) - 1, (
            'the argon sensitivity matrix has rank {0} over {1} parameters; the '
            'A/T0 identity makes {2} the only possible answer'.format(
                rank, len(names), len(names) - 1))
        null_direction = dict(zip(names, right[rank]))
        carried = {name for name, weight in null_direction.items()
                   if abs(weight) > 1.0e-6}
        assert carried == {'A', 'T0'}, (
            'the direction the anchors cannot see is over {0}, and the algebra '
            'says it is exactly {{A, T0}}'.format(sorted(carried)))

        # The direct parameter assertion is what catches it.
        assert 'T0' in ARGON_RECOMBINATION.rate_law_parameters
        reaction, entry = _selected_entry_carrying(ARGON_RECOMBINATION, traded)
        with pytest.raises(AssertionError) as raised:
            ARGON_RECOMBINATION.assert_as_shipped(reaction, entry=entry)
        message = str(raised.value)
        assert 'wrong RATE LAW PARAMETER' in message and 'T0' in message, (
            'nothing in this file reads T0:\n{0}'.format(message))
        assert 'wrong RATE MAGNITUDE' not in message, (
            'an anchor reported the T0 trade, which contradicts the algebra this '
            'test asserts:\n{0}'.format(message))

    @pytest.mark.parametrize('bound', ['Tmin', 'Tmax'])
    def test_a_non_finite_validity_bound_is_refused(self, bound):
        """``Tmin = NaN`` and ``Tmax = NaN``, each on its own.

        The window comparison used to be ``abs(found - expected) > 1e-9 *
        abs(expected)`` -- an "is it wrong" test, whose failure branch is the one
        NaN takes. So a NaN bound passed. Nothing else caught it either: measured,
        ``max(1e4, nan)`` is ``1e4`` and ``min(1e8, nan)`` is ``1e8``, so the Te
        response clamps to the unchanged band and returns its usual answer, while
        production's ``is_temperature_valid(20000)`` returns ``False`` on that
        same object -- an entry RMG will refuse to use, which this file called
        correct.

        The two bounds are checked independently because a check that only ever
        sees them together cannot tell which one it is reading.
        """
        window = dict(Tmin=(1.0e4, 'K'), Tmax=(1.0e8, 'K'))
        window[bound] = (float('nan'), 'K')
        kinetics = TwoTemperaturePlasma(
            A=(_SVS_ARGON['A_rad'], 'cm^3/(molecule*s)'), n=-_SVS_ARGON['eta'],
            Ea_g=(0.0, 'kJ/mol'), Ea_e=(0.0, 'kJ/mol'), T0=(1.0e4, 'K'),
            electrons=-1, **window)

        # The premise: production refuses this object, and the Te response does not.
        assert kinetics.is_temperature_valid(2.0e4) is False, (
            'production now accepts a NaN {0}, so this test is measuring something '
            'other than the divergence it was written for'.format(bound))
        assert electron_temperature_response(kinetics) == FALLS, (
            'the Te response now reports a NaN {0}, so the window assertion is no '
            'longer the only thing that can see it and this test should say '
            'so'.format(bound))

        reaction, entry = _selected_entry_carrying(ARGON_RECOMBINATION, kinetics)
        with pytest.raises(AssertionError) as raised:
            ARGON_RECOMBINATION.assert_as_shipped(reaction, entry=entry)
        message = str(raised.value)
        assert 'wrong TEMPERATURE WINDOW' in message and bound in message, (
            'a NaN {0} is not refused by name:\n{1}'.format(bound, message))
        other = 'Tmax' if bound == 'Tmin' else 'Tmin'
        assert '{0}=nan'.format(other) not in message, (
            'the refusal blames {0} as well, so the two bounds are not being read '
            'independently:\n{1}'.format(other, message))

    def test_the_nan_polarity_census_is_true_of_the_real_comparisons(self):
        """:data:`NAN_POLARITY_CENSUS` is a claim about behaviour, so it is driven.

        The defect repaired this round was one comparison of the wrong polarity,
        and the class is broader than the line: a comparison phrased as "is this
        value wrong" hands NaN a free pass, and one phrased as "is this value
        right" refuses it. The census records which of this file's numeric
        comparisons are which; this drives each of them with NaN and checks the
        recorded verdict, so a comparison that flips polarity turns this red
        rather than quietly joining the other class.
        """
        nan = float('nan')
        verdicts = {}

        anchor = RateAnchor(2.0e4, 1.0, 'a fixture')
        verdicts['RateAnchor.agrees'] = (
            'REFUSES' if not anchor.agrees(nan) else 'PASSES')

        verdicts['ShippedEntry.rate_law_parameters via _matches_within'] = (
            'REFUSES' if not _matches_within(nan, 1.0e4, 1e-9) else 'PASSES')
        verdicts['ShippedEntry.temperature_window via _matches_within'] = (
            'REFUSES' if not _matches_within(nan, 1.0e4, 1e-9) else 'PASSES')

        # The clamp, driven exactly as electron_temperature_response drives it.
        clamped_low = max(TE_BAND[0], nan)
        clamped_high = min(TE_BAND[1], nan)
        verdicts['electron_temperature_response band clamp'] = (
            'PASSES' if (clamped_low == TE_BAND[0] and clamped_high == TE_BAND[1])
            else 'REFUSES')

        for name, (recorded, why) in sorted(NAN_POLARITY_CENSUS.items()):
            assert verdicts[name] == recorded, (
                '{0} is recorded as {1} on NaN and measures as {2}. The census '
                'says: {3}'.format(name, recorded, verdicts[name], why))
        assert set(verdicts) == set(NAN_POLARITY_CENSUS), (
            'the census and the comparisons driven here have diverged: {0}'.format(
                set(verdicts) ^ set(NAN_POLARITY_CENSUS)))


@pytest.mark.database
class TestTheReferenceRateLawsMatchTheDatabase(_PlasmaLibraryFixture):
    """Tie :data:`ANCHORED_ENTRIES`'s fixtures back to the entries that ship.

    The adequacy checks above run against rate laws built in this file, so that
    they measure the design of the anchor set rather than the state of a
    checkout. That makes them a fourth transcription unless something compares
    them to the real thing, which is what this does.
    """

    @pytest.mark.parametrize('case', ANCHORED_ENTRIES, ids=_anchored_entry_ids)
    def test_the_reference_rate_laws_are_what_the_database_ships(self, case):
        reference, build = case
        shipped = self._select_entry(self._database(), reference).data
        built = build()
        assert type(shipped) is type(built), (
            '{0} ships a {1} and this file\'s fixture is a {2}'.format(
                reference.name, type(shipped).__name__, type(built).__name__))
        for name in rate_law_parameter_names(built):
            held = getattr(shipped, name, None)
            found = held if isinstance(held, float) else getattr(held, 'value_si', None)
            held = getattr(built, name, None)
            wanted = held if isinstance(held, float) else getattr(held, 'value_si', None)
            assert _matches_within(found, wanted, 1e-9), (
                '{0}: the shipped entry declares {1} = {2!r} (SI) and the fixture '
                'used by the adequacy checks declares {3!r}'.format(
                    reference.name, name, found, wanted))


@pytest.mark.database
class TestTheShippedCrossSectionsAreRead:
    """The real ``ElectronCollisionPlasma`` entries, swept on their own terms.

    Until this round the only ``ElectronCollisionPlasma`` this suite ever read was
    a synthetic two-point cross-section built in a fixture, and the three shipped
    libraries that hold the real ones -- ``PlasmaAir``, ``PlasmaAlkali`` and
    ``PlasmaArgon`` -- were excluded from every sweep. The argon ionisation entry
    is the single piece of data this campaign most depends on, and it sat outside
    the end-to-end discriminator entirely.

    **A measurement corrected the reason that exclusion was given.** The comment
    on ``LIBRARIES_NOT_SWEPT`` said those three libraries "are written in
    ``ElectronCollisionPlasma``". They are not, and
    :meth:`test_the_three_excluded_libraries_are_not_what_the_exclusion_claimed`
    is that retraction made into a check. What is true is narrower and is the real
    reason they cannot simply join the other sweep: they carry unimolecular,
    bimolecular and three-body entries together, so no single order holds across
    them, and the ionisation-channel Te expectation does not hold across them
    either.

    What can be swept over them, and is, is every entry of the class this suite
    could not read at all until last round, plus the weaker but universal claim
    that no shipped plasma entry leaves the order undetermined.
    """

    #: Every library in this database that carries plasma kinetics. Named here
    #: rather than discovered, so that a library added to the database and not to
    #: this list is a gap someone has to close by hand -- discovery would absorb it
    #: silently, and the exclusion this class exists to undo was exactly such a
    #: silence.
    PLASMA_LIBRARIES = (IONISATION, RECOMBINATION, 'PlasmaAir', 'PlasmaAlkali',
                        'PlasmaArgon')

    #: The cross-section entries that must be present, by provenance. The sweep's
    #: denominator is the collection it actually walks; this is the separate claim
    #: that the collection is not quietly missing the entries the sweep exists for.
    #: ``PlasmaArgon``'s is the one this campaign depends on.
    #:
    #: **Round 84 correction.** This roll named 7 of the 9 cross-sections the
    #: database ships: ``PlasmaAlkali``'s magnesium and silicon entries were
    #: missing, and removing both from the database left every sweep green -- the
    #: same defect the roll was written to prevent, one level down. The two are
    #: added here, and :meth:`test_the_roll_of_cross_sections_is_complete` is what
    #: keeps it from happening again: it fails if the database ships a
    #: cross-section this roll does not name.
    #: Fixed, committed sample count for the band sweep below -- not something
    #: discovered at run time. Measured at this count over TE_BAND=[1e4,1e6] K
    #: (log-spaced): 8 of the 9 shipped cross-sections are monotone
    #: non-decreasing; ``PlasmaAlkali:Na`` is not -- it turns over twice near the
    #: top of the band (Te=681292 K: k=5.9553e10 -> Te=825404 K: k=5.9356e10 ->
    #: Te=1e6 K: k=5.8649e10, a ~1.6% drop from the peak). All 9, including Na,
    #: stay at or above their value at the band minimum across every sample.
    #: That is why the sweep below asserts the floor property and not strict
    #: monotonicity: the data refuses the stronger claim.
    BAND_SWEEP_SAMPLES = 25

    NAMED_CROSS_SECTIONS = {
        'PlasmaArgon': ('Ar + e- => Arp + e- + e-',),
        'PlasmaAir': ('Ar + e- => Arp + e- + e-', 'He + e- => Hep + e- + e-',
                      'Ne + e- => Nep + e- + e-'),
        'PlasmaAlkali': ('Li + e- => Lip + e- + e-', 'Na + e- => Nap + e- + e-',
                         'K + e- => Kp + e- + e-', 'Mg + e- => Mgp + e- + e-',
                         'Si + e- => Sip + e- + e-'),
    }

    @classmethod
    def setup_class(cls):
        cls.libraries_path = os.path.join(settings['database.directory'],
                                          'kinetics', 'libraries')

    def _database(self):
        database = LoadedLibraries(self.libraries_path, list(self.PLASMA_LIBRARIES))
        assert set(database.libraries) == set(self.PLASMA_LIBRARIES), (
            'loaded {0}, wanted {1} -- a library named in PLASMA_LIBRARIES did not '
            'load, so every sweep below is over fewer libraries than it '
            'claims'.format(sorted(database.libraries), sorted(self.PLASMA_LIBRARIES)))
        return database

    def _cross_sections(self, database):
        """Every shipped ``ElectronCollisionPlasma`` entry, keyed by provenance.

        The denominator is derived from the collection and then checked two ways:
        it must be non-empty, and every entry named in
        :data:`NAMED_CROSS_SECTIONS` must be in it. A sweep with only the first
        check is green on any non-empty subset -- the library could lose the argon
        entry and the sweep would still report that everything it found was fine.
        """
        found = {}
        for label in self.PLASMA_LIBRARIES:
            for key, entry in database.libraries[label].entries.items():
                if isinstance(entry.data, ElectronCollisionPlasma):
                    found['{0}:{1}'.format(label, entry.label)] = entry
        assert found, ('no ElectronCollisionPlasma entry was found in any of {0}, so '
                       'this sweep is vacuous'.format(list(self.PLASMA_LIBRARIES)))
        for label, labels in self.NAMED_CROSS_SECTIONS.items():
            for entry_label in labels:
                key = '{0}:{1}'.format(label, entry_label)
                assert key in found, (
                    'the cross-section entry {0!r} is no longer in {1} as an '
                    'ElectronCollisionPlasma, so this sweep no longer covers it.\n'
                    '  swept {2} entr(ies): {3}'.format(
                        entry_label, label, len(found), sorted(found)))
        return found

    def test_the_roll_of_cross_sections_is_complete(self):
        """The roll must name every cross-section the database ships, not merely
        some of them.

        :data:`NAMED_CROSS_SECTIONS` is what stops the sweeps going green over a
        shrinking denominator, and until round 84 it named 7 of 9 -- so the two it
        did not name could be deleted from the database with every sweep still
        reporting that everything it found was fine. A roll that is allowed to be
        incomplete is not a roll.

        **What this costs, said plainly.** The sweeps stay growth-tolerant; this
        check does not. A cross-section added to the database turns this one red
        until someone names it here. That is deliberate and it is the only
        available shape: a roll derived from the database cannot detect the
        database losing an entry, because it would lose it too. One line of
        maintenance per new cross-section buys the property that no cross-section
        can leave silently.
        """
        database = self._database()
        found = set()
        for label in self.PLASMA_LIBRARIES:
            for entry in database.libraries[label].entries.values():
                if isinstance(entry.data, ElectronCollisionPlasma):
                    found.add('{0}:{1}'.format(label, entry.label))
        named = {'{0}:{1}'.format(label, entry_label)
                 for label, labels in self.NAMED_CROSS_SECTIONS.items()
                 for entry_label in labels}

        unnamed = sorted(found - named)
        assert not unnamed, (
            'the database ships {0} ElectronCollisionPlasma entr(ies) that this roll '
            'does not name:\n    {1}\nEvery sweep in this class walks them, and '
            'because the roll does not name them, DELETING them from the database '
            'would leave every sweep green. Add each to NAMED_CROSS_SECTIONS.'.format(
                len(unnamed), '\n    '.join(unnamed)))
        missing = sorted(named - found)
        assert not missing, (
            'this roll names {0} cross-section(s) the database does not ship:\n    '
            '{1}'.format(len(missing), '\n    '.join(missing)))
        assert len(found) == 9, (
            'the database ships {0} cross-sections, and 9 were measured on database '
            '96f2afa4a. The count is pinned so that a change of denominator is a '
            'decision someone makes rather than a number that drifts: {1}'.format(
                len(found), sorted(found)))

    def test_the_three_excluded_libraries_are_not_what_the_exclusion_claimed(self):
        """The retraction, as a check rather than as a corrected comment.

        The claim was that ``PlasmaAir``, ``PlasmaAlkali`` and ``PlasmaArgon`` "are
        written in ``ElectronCollisionPlasma``". Measured on database 96f2afa4a:
        only ``PlasmaArgon``, which has one entry, is. The other two are mostly
        ``Arrhenius`` and ``TwoTemperaturePlasma`` with a handful of ``ThirdBody``,
        and a few cross-sections among them.

        This matters beyond tidiness: the claim is what justified not sweeping
        them, and it made the exclusion sound like a limitation of this suite's
        readers rather than what it is -- those libraries are chemically
        heterogeneous, so a single-order or single-response sweep over them would
        be false.
        """
        database = self._database()
        census = {}
        for label in ('PlasmaAir', 'PlasmaAlkali', 'PlasmaArgon'):
            classes = [type(entry.data).__name__
                       for entry in database.libraries[label].entries.values()]
            census[label] = {name: classes.count(name) for name in sorted(set(classes))}

        assert set(census['PlasmaArgon']) == {'ElectronCollisionPlasma'}, census
        for label in ('PlasmaAir', 'PlasmaAlkali'):
            assert len(census[label]) > 1, (
                '{0} is written in a single kinetics class after all, so the claim this '
                'check retracts may have become true: {1}'.format(label, census))
            assert census[label].get('ElectronCollisionPlasma', 0) < sum(census[label].values()), (
                '{0} IS written wholly in ElectronCollisionPlasma: {1}'.format(label, census))

        orders = set()
        for label in ('PlasmaAir', 'PlasmaAlkali'):
            for entry in database.libraries[label].entries.values():
                order = rate_order(entry.data)
                if isinstance(order, int):
                    orders.add(order)
        assert len(orders) > 1, (
            'every entry in the two heterogeneous libraries reads at one order {0}, so '
            'the real reason they are swept separately -- that no single order holds '
            'across them -- no longer applies and they could join the other '
            'sweep'.format(sorted(orders)))

    def test_every_shipped_cross_section_reads_as_second_order(self):
        """The engine's own reading of a tabulated cross-section, on real data.

        ``ElectronCollisionPlasma`` states no rate-coefficient units at all -- its
        rate is a Maxwellian average over a tabulated sigma(E) -- so the order
        comes from the engine fallback, :func:`get_plasma_rate_order`. Until this
        round that fallback was exercised only against a two-point cross-section
        invented in a fixture. These are the shipped ones.
        """
        database = self._database()
        measured = {key: rate_order(entry.data)
                    for key, entry in self._cross_sections(database).items()}
        assert set(measured.values()) == {2}, measured

    def test_every_shipped_cross_section_rises_with_te(self):
        """All of them are electron-impact ionisation, ``X + e- => X+ + 2 e-``.

        That is a threshold process: below the ionisation energy the cross-section
        is zero, so the Maxwellian average grows steeply with Te. A shipped
        cross-section that fell with Te would not be an ionisation channel, and
        this is the check that would say so.
        """
        database = self._database()
        measured = {key: electron_temperature_response(entry.data)
                    for key, entry in self._cross_sections(database).items()}
        assert set(measured.values()) == {RISES}, measured

    def _te_window(self, kinetics):
        """The same window ``electron_temperature_response`` computes: TE_BAND
        intersected with the rate's own declared ``[Tmin, Tmax]``."""
        low, high = TE_BAND
        declared_low = getattr(getattr(kinetics, 'Tmin', None), 'value_si', None)
        declared_high = getattr(getattr(kinetics, 'Tmax', None), 'value_si', None)
        if declared_low is not None:
            low = max(low, declared_low)
        if declared_high is not None:
            high = min(high, declared_high)
        return low, high

    def test_every_shipped_cross_section_never_drops_below_its_band_minimum(self):
        """The real sweep the two-point endpoint ratio above cannot see, strengthened
        as far as the data actually licenses -- and no further.

        **Measured** (``BAND_SWEEP_SAMPLES`` log-spaced ``Te`` samples across
        ``TE_BAND`` intersected with each entry's declared window, evaluated through
        :func:`evaluate_as_the_solver_does` exactly as
        ``test_every_shipped_cross_section_rises_with_te`` does): 8 of the 9 shipped
        ``ElectronCollisionPlasma`` entries are monotone non-decreasing across every
        sample. ``PlasmaAlkali:Na`` is not: it turns over twice near the top of the
        band, ~1.6% down from its peak by ``Te = 1e6`` K. So a full monotone-sweep
        assertion is refused by the data -- asserting it would be tuning the check to
        pass on a measurement that says it shouldn't, which is exactly what this
        suite exists to not do.

        What the data DOES license, and what this asserts instead: every sampled
        rate, on every shipped entry including ``Na``, is at or above the rate at the
        band minimum. That is weaker than monotonicity -- it tolerates a mid-band
        dip that does not undercut the low-Te floor -- but it is real: it is derived
        from ``BAND_SWEEP_SAMPLES`` independent samples rather than the two endpoints
        the ratio check reads, so a shipped fit that dips anywhere in the interior
        and comes back up, without ever falling below its own low-Te value, is
        exactly what this catches and the endpoint ratio does not: the endpoint
        ratio would report the same RISES verdict for such a fit unchanged.

        **Blind spot, stated plainly and measured, not asserted here:** this sweep
        does not read strict monotonicity. A synthetic sequence
        ``[10, 11, 10.5, 12, 15, 20]`` -- a dip at index 2 that never goes below the
        first sample -- passes this check (``all(k >= k[0])`` is true) while failing
        a monotone-non-decreasing check. Probed directly against the check's own
        logic in ``$TMPDIR`` (not committed): the floor assertion returned green on
        that sequence and a monotone-sweep assertion applied to the same sequence
        returned red. That is the property this test does not read.
        """
        database = self._database()
        failures = {}
        for key, entry in self._cross_sections(database).items():
            low, high = self._te_window(entry.data)
            assert high > low, '{0}: no usable Te window'.format(key)
            Tes = np.logspace(np.log10(low), np.log10(high), self.BAND_SWEEP_SAMPLES)
            ks = [evaluate_as_the_solver_does(entry.data, Te) for Te in Tes]
            band_min = ks[0]
            below = [(Te, k) for Te, k in zip(Tes, ks) if k < band_min]
            if below:
                failures[key] = below
        assert not failures, (
            'these shipped cross-sections drop below their own band-minimum rate '
            'somewhere in the sweep:\n  {0}'.format(
                '\n  '.join('{0}: {1}'.format(k, v) for k, v in failures.items())))

    def test_no_shipped_plasma_entry_leaves_the_order_undetermined(self):
        """The weaker claim, over every entry of every plasma library.

        "Undetermined" is the outcome that means no check happened, so a shipped
        entry that produces one is a hole in every order-based check downstream.
        This sweep is what turned the composite-kinetics gap from a code reading
        into a measurement: before ``_rate_coefficient_units`` learned where
        ``ThirdBody`` keeps its units, 16 shipped entries failed here.
        """
        self._assert_no_entry_leaves_the_order_undetermined(self._database())

    def _assert_no_entry_leaves_the_order_undetermined(self, database):
        """The sweep itself, factored out so a fabricated entry can be put through
        **this** code and not through a second copy of it.

        A demonstration that a mutated entry "would be refused" is worth nothing if
        it re-implements the refusal: the mirror shape this campaign keeps finding
        is N normative copies compared to each other and none to the implementation.
        :meth:`test_a_mixed_composite_entry_is_refused_by_this_same_sweep` calls
        exactly this method.
        """
        undetermined = {}
        swept = 0
        for label in self.PLASMA_LIBRARIES:
            entries = list(database.libraries[label].entries.values())
            assert entries, '{0} loaded no entries, so this sweep is vacuous on it'.format(label)
            for entry in entries:
                swept += 1
                order = rate_order(entry.data)
                if not isinstance(order, int):
                    undetermined['{0}:{1}'.format(label, entry.label)] = str(order)
        assert swept == sum(len(database.libraries[label].entries)
                            for label in self.PLASMA_LIBRARIES), 'the sweep skipped entries'
        assert not undetermined, (
            '{0} of {1} shipped plasma entries state no order, so every order-based '
            'check is silently switched off for them:\n  {2}'.format(
                len(undetermined), swept,
                '\n  '.join('{0}  {1}'.format(k, v) for k, v in sorted(undetermined.items()))))

    def test_a_mixed_composite_entry_is_refused_by_this_same_sweep(self):
        """The consequence of the composite repair, on the real sweep and real data.

        The reader-level check
        (``test_a_composite_does_not_inherit_the_order_of_the_half_it_could_read``)
        says what :func:`rate_order` now returns. This says what that buys: an entry
        whose rate law is a composite with one unreadable component is REFUSED by
        the shipped-entry sweep, where before the repair it inherited the readable
        component's order and passed.

        The entry is fabricated over a real loaded library rather than authored into
        the database, because the database is read-only to this branch -- and the
        fabrication is the nearest wrong thing to the accepted case, not an
        obviously broken one: a ``MultiArrhenius`` whose first component is a
        perfectly good second-order ``Arrhenius``. Before the repair this sweep was
        green on it.

        The restore is asserted, not assumed: this method mutates a live loaded
        entry, and a sweep left looking at a fabricated rate law would poison every
        later check in the class.
        """
        database = self._database()
        label = self.PLASMA_LIBRARIES[0]
        victim = sorted(database.libraries[label].entries)[0]
        entry = database.libraries[label].entries[victim]
        original = entry.data

        readable = Arrhenius(A=(1.0e10, 'm^3/(mol*s)'), n=0.0, Ea=(0.0, 'kJ/mol'))
        mixed = MultiArrhenius(arrhenius=[readable, Arrhenius()])
        assert rate_order(readable) == 2, 'the control is not a control'

        try:
            entry.data = mixed
            with pytest.raises(AssertionError) as raised:
                self._assert_no_entry_leaves_the_order_undetermined(database)
            message = str(raised.value)
            assert COMPONENT_WITHOUT_UNITS in message, (
                'the sweep refused the entry but not for this reason:\n{0}'.format(message))
            assert '{0}:{1}'.format(label, entry.label) in message, (
                'the refusal does not name the entry it is about:\n{0}'.format(message))
        finally:
            entry.data = original
        assert entry.data is original, 'the fabricated rate law was left on the entry'
        self._assert_no_entry_leaves_the_order_undetermined(database)


@pytest.mark.database
class TestTheLibraryLoadLandmineIsNotReachedFromHere(_PlasmaLibraryFixture):
    """A landmine adjacent to this repair, measured and NOT fixed by it.

    ``KineticsLibrary.check_for_duplicates`` runs inside ``load`` and raises
    ``DatabaseError`` on any unmarked duplicate, deciding duplication with
    ``Reaction.is_isomorphic``. So it looked as though repairing that predicate
    would also let one library carry both channels. Measured, it does not, and the
    reason is upstream of anything this change touches: the check runs over
    ``entry.item``, which ``KineticsLibrary.load_entry`` builds as a plain
    :class:`Reaction` **with no owner recorded at all**. No owner means no
    placement declaration, which means the net-derived rule -- under which the two
    channels are ``(0, 1)`` and ``(1, 0)``, exact mirrors, and therefore still
    duplicates.

    The owner only appears later, in ``get_library_reactions``, which wraps each
    entry in a :class:`LibraryReaction` whose ``family`` is the library label.
    Everything downstream of that point -- the model builder, the reactor, both
    writers -- sees the declaration and is repaired. The library loader is the one
    consumer that sees the reaction before its owner is attached.

    The verdict is also **unchanged** by this commit rather than broken by it: the
    old code compared ``electrons == -electrons`` (``+1`` against ``-(-1)``, true)
    and the new code compares ``(0, 1)`` against the mirror of ``(1, 0)`` (equal,
    true). Same answer, same reason.

    Fixing it means giving library entries their owner at load, and it needs a
    design decision this ticket has no mandate for: the placement registry is
    keyed per owner, so one library cannot declare two different placements even
    once the owner is recorded. That is a separate ticket, and the strict xfail
    below is what will tell whoever takes it that they have succeeded.

    The library here is not a fixture: its entries are the shipped
    ``PlasmaElectronImpactIonization`` and ``PlasmaRadiativeRecombination``
    entries, loaded out of the database by RMG's own loader.
    """

    def _merged_library(self, database):
        """One library carrying both lithium channels, and nothing else.

        Built by *selecting* the two entries it means, one out of each shipped
        library. It used to sweep in every entry of both and then assert the
        total was two -- which was a statement about the database's size, and
        stopped being true the moment the recombination library grew its argon
        entry. The count below is not a return of that defect: ``merged.entries``
        is a collection this helper filled two lines earlier from two named
        selections, so counting it asserts that the two selections did not land
        on the same entry, which is a property of this helper alone.
        """
        merged = KineticsLibrary(label='PlasmaBothChannels')
        merged.entries = {}
        for index, label in enumerate((IONISATION, RECOMBINATION), start=1):
            entry = self._select_entry(database, LITHIUM_CHANNEL[label])
            entry.index = index
            merged.entries[index] = entry
        assert len(merged.entries) == 2, 'the two channel selections collided on one entry'
        return merged

    def test_a_loaded_entry_carries_no_owner(self):
        """The measurement the whole limitation rests on.

        The denominator is asserted, and asserted the way the rest of this file
        now does it: derived from the collection walked, never a frozen count.
        Without it this check was an unguarded ``for`` over a collection it did
        not own -- green on a library that loaded nothing at all, which is the
        state in which the limitation it claims to measure is not measured at all.

        Two separate claims, because they fail separately. *Coverage*: every
        library walked yielded entries, and the number of assertions made equals
        the number of entries there were to make them about. *Named presence*: the
        two channels this file is about are among them, found by provenance -- a
        sweep of a library that had lost the lithium entries but kept the argon
        one would satisfy coverage and measure the wrong thing.
        """
        database = self._database()
        seen = 0
        for label in (IONISATION, RECOMBINATION):
            entries = list(database.libraries[label].entries.values())
            assert entries, (
                '{0} loaded no entries, so this measurement was not made on it'.format(label))
            for entry in entries:
                assert getattr(entry.item, 'family', None) is None
                assert get_placement_declaration(entry.item) is None
                seen += 1
        assert seen == sum(len(database.libraries[label].entries)
                           for label in (IONISATION, RECOMBINATION)), (
            'the sweep made fewer measurements than there were entries to measure')
        for label in (IONISATION, RECOMBINATION):
            named = self._select_entry(database, LITHIUM_CHANNEL[label])
            assert getattr(named.item, 'family', None) is None, (
                'the {0} entry this file is about now carries an owner at load, which '
                'is the ticket the xfail below is waiting for'.format(label))

    def test_a_loaded_entry_carries_no_kinetics_either(self):
        """The owner is not the only thing ``entry.item`` is missing.

        The loader parks the rate on ``entry.data`` and only
        ``get_library_reactions`` puts item and data together, so at this stage a
        candidate can state neither its owner nor its channel. Measured rather
        than assumed, because ``_select_entry`` is built on it: reading
        ``reaction.kinetics`` over raw entries would have made the rate-order
        discriminator vacuously ``None == None``-free -- it would match nothing,
        loudly, rather than match everything, quietly, but either way the reason
        belongs in a check and not in a comment.
        """
        database = self._database()
        seen = 0
        for label in (IONISATION, RECOMBINATION):
            for entry in database.libraries[label].entries.values():
                assert entry.item.kinetics is None
                assert entry.data is not None
                assert rate_order(entry.data) == 2
                seen += 1
        assert seen, 'neither library yielded an entry, so nothing above was measured'

    def test_the_owner_appears_when_the_entry_becomes_a_library_reaction(self):
        """And from there on, the placement is available and this repair applies."""
        database = self._database()
        placements = {}
        for label in (IONISATION, RECOMBINATION):
            reaction = self._select(database, LITHIUM_CHANNEL[label])
            assert reaction.family == label
            placements[label] = get_electron_placement_counts(reaction)
        assert placements == {IONISATION: (1, 2), RECOMBINATION: (1, 0)}

    @pytest.mark.xfail(strict=True, raises=DatabaseError, reason=(
        'Not reached by this repair: check_for_duplicates runs over entry.item, which '
        'carries no owner, so both channels fall back to the net rule and stay mirrors. '
        'Fixing it requires recording the owner at load, and a decision about a registry '
        'that is keyed one placement per owner. When this starts passing, that ticket has '
        'landed and this marker should be removed. `raises=DatabaseError` pins WHICH '
        'failure is the expected one: without it the marker was satisfied by an '
        'AssertionError raised in the set-up helper two statements before '
        'check_for_duplicates was ever called, so the check reported xfailed -- green '
        'to any reader -- while measuring nothing at all, and its success signal, the '
        'strict xfail flipping to a pass when that other ticket lands, was silently '
        'disabled.'))
    def test_one_library_carrying_both_channels_can_be_loaded(self):
        self._merged_library(self._database()).check_for_duplicates()
