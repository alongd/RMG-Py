"""Recertified closed windows, chain-end contexts and the short-chain ladder
for the declared v1 PS domain (I-137).

For one record member (a channel member: one mapped outcome, with every RMG
mapping of it) on a boundary witness (``rmgpy.kmc.boundaries``), this module
builds a CLOSED WINDOW that provably contains everything RMG read to select
the member's rate, stable on every extension of the witness, or refuses by
name.  It enumerates the chain-end contexts of the member and the finite
short-chain ladder per family.  No templates or weights are cut here.

Support (item 1)
----------------
The support of a member is the union, over every mapping of the member (A+A
template slot orders and reactant resonance variants included, none
skipped), of

* the reacting atoms: labelled atoms (H-only centres seed with their labelled
  H), every atom of ``Move.mapped`` and every atom of the recipe's edits;
* for every predicate RMG tried in every descent (root checks, matched AND
  failed children, logic components in RMG's evaluation order): its anchors,
  the candidates of its free labels (when bounded, below), its ball (the
  predicate's radius around anchors and candidates, ``certificates``), the
  images of its embeddings and the labelled atoms RMG ignored;
* for a reverse record nothing else: the certificate is built on the family
  template side (the record's products), which shares the witness's atom-id
  namespace, so the support is transported by the identity map (``identity``
  is checked: every support atom is a witness atom).

Radius unit: ordinary bonds as in the certificate (``RADIUS_UNIT``); explicit
hydrogens are graph atoms, so a ball reaches the H of its outer atoms only
when the radius does; closure adds every H of an included heavy atom.  The
boundary signature (full bond lists of the ball rows) names the bonds that
leave a ball: in the declared domain (C and H only, no charge or lone pair)
RMG's atom type of a ball atom is a function of its own electrons and bond
orders, never of the far atom's element, so the far atom is not read; the
bond's ORDER is, and the closure makes it a window property (any bond of
order > 1 lies in a pi system or a ring the closure includes whole; checked:
``boundary_bond_problems``).

Not atom sets (each with an explicit disposition, ``Support.dispositions``):

* free-label searches (R_Recombination's two ``*``): RMG scans the whole
  structure.  ``bounded_by_domain`` when no atom an extension can add (an
  interior repeat unit of the provider, its atom types computed by RMG) is a
  specific case of the free group atom (so no extension adds a candidate:
  structural, not observed) AND every candidate is a reacting atom; a
  candidate outside the reacting atoms (a remote radical site) makes the
  selection depend on the whole species: ``whole_species`` -> refused; an
  extension able to add candidates -> refused (``unbounded_search``);
* group multiplicity and radical-count guards (whole molecule): a context
  condition, the species' unpaired-electron count under the declared
  maximum-multiplicity convention; extension invariant because every provider
  fragment is closed shell (checked per fragment); varied only by the radical
  inventory of the contexts, which ``context_assignments`` always enumerates;
  the guard values are part of the certificate view;
* element-count guards: monotone under extension (the net element delta of a
  provider step is non-negative, checked): a satisfied minimum stays
  satisfied; a failed one with zero delta stays failed; a failed one with a
  positive delta could flip at some depth -> refused;
* ``inRing``: closure includes every ring <= ``RING_LIMIT`` through a window
  atom;
* the depository search (whole species): a depository match makes the
  certificate ``whole_species`` -> refused (no window); a rule estimate needs
  a non-refused depository exclusion;
* the multiplicity convention: declared, checked on every species by
  ``domain.check``;
* the channel kind (``kept`` | ``synthetic:unchanged_species``) and the
  identical-reactant pairing: whole-species identity verdicts (an extension of
  one participant changes them), not reads of the selection; both kinds are
  rated by the family's estimate for the template, so the kind is reported
  as a disposition and left out of the certificate view.

Closure (item 2)
----------------
On BOTH reaction sides (the window is the union): every ring <= ``RING_LIMIT``
through a window atom, whole (a reacting atom in a larger ring -> refused; a
larger ring touching the window without a reacting atom is refused too: the
declared domain admits six-membered phenyl rings only, so such a ring is a
domain error, never silently truncated); the whole conjugated pi / radical
system (atoms with an unpaired electron or a bond of order > 1, joined by
bonds) touching the window; every H of an included heavy atom; repeated to a
fixed point.

Recertified fixed point (item 3)
--------------------------------
``fixed_point``: observe (certify, support, close); while an artificial port
(port atom or a declared cap) lies in the window, extend every such port by
one repeat (``boundaries.extend``), recertify and re-close.  Then the checks
of the round (``check_groups``, every round, no option turns one off): EVERY
artificial port alone at +1 and at +2 repeats (``CHECK_STEPS``), and all
artificial ports jointly at +1 and at +2; each recertified and closed.  The
fixed point needs every check to leave the window's marked graph (identity
correspondence: surviving atoms keep their ids) and the certificate view
(``certificate_view``) unchanged.  A failed independent check extends its
port by one; a failure seen only jointly extends every port by one; the
round repeats.  Pins: ports in vector order
(participant, port id); ``PORT_DEPTH_LIMIT`` repeats per port; at most
``FIXED_POINT_RETRIES`` re-extensions after failed checks; the heavy-atom
ceiling ``PROXY_PADDING_MAX_HEAVY_ATOMS``.  Statuses: ``fixed_point``,
``no_fixed_point`` (a port is still in the window at its limit, or windows
keep changing past the retries), ``ceiling_exceeded`` (never truncated),
``certificate_unstable`` (the certificate keeps changing past the retries),
and ``refused:<reason>`` from the observation.

Span bound (item 4)
-------------------
Per family, from the declared domain and the family's tree constants
(``family_constants``): ``r`` the largest predicate radius of the tree,
``D`` the largest bond distance between labelled heavy atoms of one
connected template structure (reactants, and products after the recipe).
Every support atom lies within ``r`` bonds of a reacting atom, so the
support's backbone extent is at most ``D + 2r`` bonds.  Closure extends the
backbone by at most ``CLOSURE_EXTENSION = 1`` bond per side: rings are
pendant phenyls (any ring atom's backbone carbon is nearer the centre than
the ring atom), radical pi systems are a backbone carbon plus its pendant ring
(allylic radicals are outside the domain), and the only backbone pi system
is an end alkene (two carbons).  Hence a window spans at most
``A = D + 2(r + 1) + 1`` backbone atoms of one chain, i.e. at most
``ceil((A + 1) / 2)`` PS repeat units (two backbone carbons per repeat, any
phase); the bound on units is in each provider's own units
(``repeats_bound``: a provider with ``k`` backbone atoms per repeat,
``ceil((A + k - 1) / k)``; the alkane's one-carbon repeat: ``A``).  A
window exceeding its bound fails the run (``SpanBoundError``).

Ladder (item 6)
---------------
A chain of DP ``n`` has ``2n - h + t`` backbone carbons (``h``/``t``: the
head/tail starts/ends on the other repeat phase, ``PsChainSpec.head_shift``
/``tail_shift``), at least ``2n - 1``.  A window can hold both physical ends
of a chain iff the chain's backbone fits its span: ``L = floor((A + 1) / 2)``
(the largest n with 2n - 1 <= A; in-phase chains need 2n <= A, so the
repeat phase and the end groups decide which chains of DP L reach).  DP
``L + 1`` (``2L + 1 > A`` backbone carbons in every phase) is the
non-overlap control.  DP 1 includes the one-carbon species (toluene, benzyl
radical, methane, methyl radical).  Bimolecular directions: the participant
kinds (radical, closed shell, alkene) each direction needs come from the
family's own templates (``pair_roster``); every DP pair (a, b) with
1 <= a, b <= L + 1 is planned for every kind combination in its direction
(``ladder_plan``) and coverage is checked per direction
(``ladder_direction_gaps``).
This DP grid runs on one representative per participant kind
(``kind_specs``); the generated roster below adds every reacting-site
environment class.

Generated roster (reacting-site environment classes)
----------------------------------------------------
A reacting-site environment class is a reacting site (a site kind: every
radical site kind, every closed-shell carbon, every alkene end carbon) on a
participant pattern (the declared features the participant carries: its two
end groups and its interior radicals; a one-carbon or small species), in one
template slot of one direction, with the smallest DP at which RMG generates
an in-domain member carrying it.  Why the classes exhaust ``domain.py``'s
declared features:

* the species ``domain.check`` accepts are linear PS chains (a simple path
  of backbone carbons, phenyl pendants, every other substituent H), the
  one-carbon chains and the small species ``SMALL_SPECIES``;
* a chain's features are its two termini, each one of ``END_GROUPS`` (the
  only declared end groups; unsaturation is allowed at an end only), and its
  interior radical sites, each of a ``RADICAL_KINDS`` kind that can sit on an
  interior backbone carbon (``INTERIOR_RADICALS``: benzylic, secondary; a
  primary site is a terminus, a methyl site a one-carbon chain, an H atom a
  small species), at most ``MAX_RADICAL_SITES`` radical sites in all; no
  other feature is declared (``declared_site_kinds`` maps every declared key
  to its site kinds and raises on an unmapped one);
* so ``participant_patterns`` (every unordered pair of end groups, every
  multiset of interior radical kinds within the radical budget, the
  one-carbon and small species) is every declared participant up to the
  placement of its features along the chain;
* a reacting site is one carbon of such a participant (or the H atom / H2),
  typed by ``site_kind`` from the same reading as ``domain`` (radical,
  terminus or interior, phenyl, double bond); the features that can lie
  within the span bound of it are among the pattern's, so (site kind,
  pattern) is the class, and the site's template slot (label set) is read
  from RMG's own mappings (``probe_pattern``);
* placement: ``pattern_specs`` puts the interior radicals as near the head
  and as near the tail as the domain allows (the canonical placements), and
  the probe runs each pattern from its smallest constructible in-domain DP
  through ``probe_extra`` more repeats (every distance constraint of the
  domain is a lower bound of at most ``RADICAL_SEPARATION + 1`` bonds between
  features and the event's atoms, which span at most ``D``), at most L + 1;
  the class's smallest DP is the first at which an in-domain member carries
  it (a class can need more repeats than its species: the product radical
  of an H abstraction must sit more than ``RADICAL_SEPARATION`` bonds from
  the species' other radical).

What the class product does not vary: the distance of each feature from the
site beyond the smallest DP.  Those distances are the DP grid's (every DP
pair up to L + 1 on the kind representatives) and the contexts' (every
declared end at every near port of the long artificial chains); a feature
beyond the span bound acts only through whole-species reads (the radical
inventory), which the far-port radical ends carry.

The probe pairs a pattern, in a bimolecular direction, with the cheapest
declared species that fills the other slot in an in-domain member
(``probe_partners``): RMG matches a template slot to one reactant whatever
fills the other slot, and the partner's own products do not depend on the
pattern.  Ladder items (``class_items``): per unimolecular direction every
class representative alone; per bimolecular direction the product of the
slot-0 and slot-1 class representatives, each at its smallest DP; every class
must reach the plan (``ladder_class_gaps``).  ``ladder_plan`` raises when a
family's class roster or (bimolecular) pair roster was not computed.

Class equivalence and the reduced class plan (rework 3)
-------------------------------------------------------
The class product duplicates rate-selection environments; the plan expands
distinct ones only, and each reduction is a checked claim.

* Closing (``close_classes``): each class is closed ONCE against its
  closing partner (``closing_partners``: per slot the cheapest declared
  species filling the other slot whose member reaches a fixed point; the
  probe partner H / H2 ends ``refused:guard: element_count``): every member
  whose slot carries the class's site on the class's spec runs through
  ``run_member``.  The closing records the members' results, pathways, own
  leaves and statuses.
* Equivalence (``classes_merge``, ``merge_classes``): two classes of one
  direction and slot merge iff they were closed against the same partner,
  have the same pathways (leaves, checked explicitly), the same member
  statuses, at least one fixed point, and every fixed-point member of each
  has a member of the other whose closed window is the same marked graph
  (``class_window_isomorphisms``: id-free reacting marks, the windows sit on
  different witnesses) with an equal certificate view under that
  correspondence (reacting atoms transported).  Windows are bounded by the
  span bound, so classes differing only beyond it merge; a different leaf
  never merges.  Each class maps to the cheapest class it equals (smallest
  DP, fewest heavy atoms): the representative (``ClassMerge``).
* Merged-root trees (Disproportionation, R_Recombination: one
  auto-generated ``Root`` whose predicates read both reactants): the
  closing partner's descent decides which predicates read the class's
  participant, so the window seen with one partner does not bound what
  another partner's descent reads.  Every predicate of any descent reads
  within ``r`` bonds of the labelled atoms, so these directions also
  require the windows widened to ``r`` heavy bonds around the participant's
  reacting atoms and closed (``ball_window``) to be equal, and equal
  whole-species
  reads (unpaired electrons; element counts capped at the family's largest
  guard minimum).  They keep the PAIR PRODUCT over the representatives.
* Per-reactant trees (``per_reactant_tree``, read from the loaded family:
  not auto-generated, every forward-template reactant group its own
  top-level node; H_Abstraction and R_Addition_MultipleBond forward): RMG
  descends each reactant's tree on the reactant that carries its labels
  alone (``get_reaction_template`` matches each template group to one
  reactant), so a participant's descent, support, window and leaf do not
  depend on its partner; the rate is then looked up by the tuple of leaves
  (rule, averaging, parent walk: functions of the tuple).  Whole-species
  reads of a group (multiplicity, radical and element counts) read the
  group's own reactant.  Two reads see both species: the depository match
  (a whole-reaction isomorphism) and the channel kind (whole-species
  identity, a disposition outside the view).  Hence the cover: every class
  is closed (its representative's closing reached a fixed point, or it is
  expanded unreduced: a class whose closing reached none, every member
  refused by status, pairs with every representative of the other slot);
  one pair of representative specs per leaf tuple that occurs (slot-0
  leaves x slot-1 leaves, ``leaf_cover``; each pair item names its tuples
  and ``leaf_tuple_gaps`` needs a fixed-point member on each); and the
  product of the classes whose species is a training-reaction species
  (``training_keys``: the only pairs a depository match can read).
  Unimolecular directions have no pair product and are not closed: each
  class's own run (every representative spec alone) is its one closing.
* Expansion over representatives (``site_select``, ``LadderItem.sites``): a
  class pair expands only the members carrying one of its specs'
  representative sites in each slot, and a unimolecular class run only the
  members carrying one of its spec's own classes; every other member of the
  same species carries a class run on its own (smallest-DP) spec or merged
  into a representative.  A class item on the species pair of a grid item
  runs every member.
* Pre-filter (``pair_prefilter``): a pair is dropped before expansion only
  when a one-product template (``joined_templates``) puts more than
  ``MAX_RADICAL_SITES`` unpaired electrons on its single product in every
  mapping (both reactants' electrons plus the recipe's delta): ``domain``
  refuses every member.  Two-product templates split the electrons per
  participant, where the closings already put each class in domain.

Contexts (item 5)
-----------------
Per artificial port: artificial, or physical with a declared end chemistry
built by the provider (``boundaries.reconstruct_end``: caps deleted, fresh
end atoms, never a relabelled cap); the event is rematched on the context
by ``rmg_moves`` (same labels, same mapped outcome), recertified and closed.
The assignments (``context_assignments``) are the full product over all
artificial ports of all participants, whatever the certificate reads.  A
port farther than ``reach + 1`` heavy bonds from every reacting atom
(``member_near_ports``) cannot lie in any context's window (every window
atom is within ``reach`` bonds of a reacting atom, item 4), so it can act
only through whole-species reads: the radical inventory (its radical ends
stay in the product) and the element-count guards.  Its non-radical ends are
reduced to one single each, and that reduction is a checked claim: every
such single must merge with the all-artificial context (equal window, equal
certificate view: ``far_port_verdicts``), otherwise the port is promoted to
near and the product rebuilt (``far_port_promotions``).  Single ends merging
does not imply that their combinations merge, so the check is also joint:
every combination of non-radical ends over all far ports of one participant
and over all far ports (``far_port_groups``) must merge too, else those ports
are promoted.  Two contexts merge
only if their closed windows are equal as marked graphs (both sides,
participant roles, reacting atoms, physical ends, oriented exterior bonds)
and their certificate views are equal under that correspondence.

Reacting ports (policy a)
-------------------------
A member whose reacting atom sits on an artificial port is refused
(``refused:reacting_port``): the cut point is not a real event.  It stands for
one distinct mapped event per assignment of each reacting port to the
interior (the chain goes on there) or to a declared end chemistry (the chain
ends there), and each event gets its own witness
(``reacting_port_covers``): the specs with the interior ports one repeat
longer and the end ports built physical by the provider, the member
rematched by its labelled atoms' repeat coordinates (``event_signatures``).
An event is covered only if every matching member reaches a fixed point or
is shown outside the declared domain (``domain.check`` refuses one of its
species, reason stated: first the event's own reactants, then every
species of each matching member), or if RMG's family generates no mapping on
the event's atoms of its own witness under any labels (``no_event``: that
species has no such reaction, so no rate is selected for it); certificate
failures, ceiling hits, incomplete descents, port or end refusals and a
rematch that finds the atoms only under other labels are gaps
(``cover_status_verdict``).  The member is covered iff every event is
(``cover_verdict``: one event never stands in for another); an uncovered one
is a gap (``reacting_port_gaps``).
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import combinations, combinations_with_replacement, permutations, product
from typing import Any, Callable, Optional, Sequence

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import declared_chain as dc
from rmgpy.kmc import domain as dom

RING_LIMIT = 8
PORT_DEPTH_LIMIT = 6
FIXED_POINT_RETRIES = 3
CHECK_STEPS = (1, 2)
CLOSURE_EXTENSION = 1
ISOMORPHISM_LIMIT = 256
LOOP_STATUSES = ("fixed_point", "no_fixed_point", "ceiling_exceeded", "certificate_unstable")
DISPOSITIONS = ("bounded_by_domain", "context_condition", "monotone_satisfied",
                "invariant_failed", "closure", "declared", "whole_species",
                "whole_species_identity", "refused")


class WindowError(ValueError):
    """Bad input."""


class SpanBoundError(AssertionError):
    """A window exceeds its family's derived span bound (a bug or a domain
    error): fails the run."""


# ---------------------------------------------------------------------------
# provider objects (I-165): every chemistry the window machinery reads comes
# from the provider of a witness participant, never from a provider name
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class DeclaredSpec:
    """A witness spec of a declared-chain provider: the registered provider
    name and the chain (``declared_chain.ChainSpec``)."""

    provider: str
    chain: Any


class SmallSpeciesProvider(b.BoundaryProvider):
    """The small declared species (``small_build``): no ports, no repeat
    coordinates, no backbone (no span check), no interior atoms, no
    extension fragment and no growth unit."""

    name = "small"
    SPECIES = ("H", "H2")
    END_CHEMISTRIES: dict = {}

    def build(self, name):
        return small_build(name)

    def fragment(self, port, coordinate):
        raise b.ExtensionRefused("a small species has no ports")

    @classmethod
    def normalize_spec(cls, spec):
        """("small", name) of a declared small species; None for another
        shape; ValueError for a malformed one."""
        if not isinstance(spec, tuple) or not spec or spec[0] != cls.name:
            return None
        if len(spec) != 2 or spec[1] not in cls.SPECIES:
            raise ValueError(f"undeclared small species spec {spec!r} (declared: {cls.SPECIES})")
        return spec


SMALL_PROVIDER = SmallSpeciesProvider()
# the owners of the named witness-spec shapes, asked in order
SPEC_OWNERS = b.LEGACY_PROVIDERS + (SMALL_PROVIDER,)


def witness_spec(spec):
    """The normalized witness spec: a ``DeclaredSpec`` of a registered
    declared-chain provider, or a named shape normalized by its owner
    (``SPEC_OWNERS``).  Anything else, a bare ``ChainSpec`` and an
    interior-radical 5-tuple included, is refused (never built closed)."""
    if isinstance(spec, DeclaredSpec):
        provider = b.PROVIDERS.get(spec.provider)
        if not isinstance(provider, dc.DeclaredChainProvider):
            raise WindowError(f"refused:legacy_spec: {spec.provider!r} is not a registered "
                              "declared-chain provider")
        if not isinstance(spec.chain, dc.ChainSpec):
            raise WindowError(f"refused:legacy_spec: a DeclaredSpec needs a ChainSpec, got "
                              f"{type(spec.chain).__name__}")
        return spec
    for owner in SPEC_OWNERS:
        try:
            normalized = owner.normalize_spec(spec)
        except ValueError as error:
            raise WindowError(f"refused:legacy_spec: {error}") from error
        if normalized is not None:
            return normalized
    raise WindowError(f"refused:legacy_spec: {spec!r} is not a witness spec (a DeclaredSpec "
                      "or a named provider shape)")


def provider_object(provider, providers=None):
    """The provider object of a participant's provider name (or the object
    itself), looked up by name: never an enumeration of the registry."""
    if isinstance(provider, b.BoundaryProvider):
        return provider
    if provider == SMALL_PROVIDER.name:
        return SMALL_PROVIDER
    table = b.PROVIDERS if providers is None else providers
    found = table.get(provider)
    if found is None:
        raise WindowError(f"no provider {provider!r} is registered")
    return found


def _is_declared(provider) -> bool:
    return isinstance(provider, dc.DeclaredChainProvider)


def observation_providers(witness, providers=None) -> tuple:
    """The provider objects whose interior atoms and extension fragments a
    member observation reads: the distinct providers of the witness's
    participants when one of them is declared; otherwise the named pair
    (``boundaries.LEGACY_PROVIDERS``), exactly the base classifier and
    fragment facts."""
    found = []
    for participant in witness.participants:
        provider = provider_object(participant.metadata.provider, providers)
        if provider not in found:
            found.append(provider)
    if not any(_is_declared(p) for p in found):
        return b.LEGACY_PROVIDERS
    return tuple(found)


def growth_unit(witness, providers=None) -> Optional[tuple]:
    """The depository growth unit of a witness: its one chain provider's
    ``unit_formula`` (correction C-1: legacy witnesses included, PS C8H8 as
    before, the linear alkane CH2); None for a witness without a chain
    provider (the certificate default).  Two different chain providers are
    refused on a declared witness; on a legacy one (PS with alkane) the
    certificate default is kept, as before."""
    found = [provider_object(p.metadata.provider, providers) for p in witness.participants]
    chains = sorted({p.name for p in found if p.unit_formula() is not None})
    if not chains:
        return None
    if len(chains) != 1:
        if not any(_is_declared(p) for p in found):
            return None
        raise WindowError(f"refused:growth_unit: two chain providers in one witness: {chains}")
    return next(p for p in found if p.name == chains[0]).unit_formula()


# ---------------------------------------------------------------------------
# chemistry policy (I-165 S2b): the domain, certificate eligibility,
# recertification and the ring-edit generation filter of a witness, named
# ---------------------------------------------------------------------------
RING_EDITS = ("exclude", "admit")


@dataclass(frozen=True)
class PolicyDomain:
    """A species domain: ``check(molecule) -> domain.DomainVerdict`` and
    ``spin_verdict(molecule) -> domain.SpinVerdict``, under ``name``."""

    name: str
    check: Callable = field(compare=False)
    spin_verdict: Callable = field(compare=False)


@dataclass(frozen=True)
class ChemistryPolicy:
    """The chemistry a witness is evaluated under: its species ``domain``
    (``species_domain_refusal``, the event builds of a reacting-port cover,
    certificate eligibility and recertification) and ``ring_edits``
    (``"exclude"``: ``member_cases`` drops every mapping whose connectivity
    change touches a ring atom, ``ring_change_select``; ``"admit"``: none)."""

    name: str
    domain: PolicyDomain
    ring_edits: str

    def __post_init__(self):
        if self.ring_edits not in RING_EDITS:
            raise WindowError(f"ring_edits {self.ring_edits!r} is not one of {RING_EDITS}")


# the legacy chemistry exactly: the v1 PS domain (``domain.check``) and no
# ring edits (late-bound, so the module functions stay the ones called)
LEGACY_V1 = ChemistryPolicy("legacy_v1", PolicyDomain(
    "legacy_v1", lambda molecule: dom.check(molecule),
    lambda molecule: dom.spin_verdict(molecule)), "exclude")
# a declared chain under the INTERIM domain (``domain.declared_interim_check``)
# and no ring edits; a later slice replaces the domain.  This object admits C
# and H only; ``declared_interim`` binds a witness's declared elements
DECLARED_INTERIM = ChemistryPolicy("declared_interim_v1", PolicyDomain(
    "declared_interim_v1", lambda molecule: dom.declared_interim_check(molecule),
    lambda molecule: dom.spin_verdict(molecule)), "exclude")
_INTERIM_POLICIES = {}


def declared_interim(elements) -> ChemistryPolicy:
    """``DECLARED_INTERIM`` admitting ``elements`` (C and H always): the
    same spin verdict and ring edits; the domain check reads only the given
    set, and the set is part of the policy's identity: the name is
    ``declared_interim_v1[C,H,O]`` (sorted), so certificates
    (``domain_policy``, the recertification gate) and every cache keyed by
    the policy name separate element sets.  ``DECLARED_INTERIM`` itself
    (``declared_interim_v1``) for C and H; one object per name."""
    elements = tuple(sorted(set(elements) | {"C", "H"}))
    if elements == ("C", "H"):
        return DECLARED_INTERIM
    name = f"{DECLARED_INTERIM.name}[{','.join(elements)}]"
    if name not in _INTERIM_POLICIES:
        _INTERIM_POLICIES[name] = ChemistryPolicy(name, PolicyDomain(
            name, lambda molecule: dom.declared_interim_check(molecule, elements),
            DECLARED_INTERIM.domain.spin_verdict), DECLARED_INTERIM.ring_edits)
    return _INTERIM_POLICIES[name]


def policy_for(witness, providers=None) -> ChemistryPolicy:
    """``LEGACY_V1`` when every participant comes from a legacy provider (the
    named chains and the small species), ``DECLARED_INTERIM`` when a chain
    participant is declared, admitting the elements of the witness's own
    declarations (``declared_interim``; the legacy small species are H and
    H2); a witness mixing legacy and declared chain providers is refused."""
    found = [provider_object(p.metadata.provider, providers) for p in witness.participants]
    declared = [p.name for p in found if _is_declared(p)]
    legacy_names = {x.name for x in b.LEGACY_PROVIDERS}
    legacy = [p.name for p in found if p.name in legacy_names]
    if declared and legacy:
        raise WindowError(f"refused:policy: a witness mixes legacy {sorted(set(legacy))} and "
                          f"declared {sorted(set(declared))} chain providers")
    if not declared:
        return LEGACY_V1
    return declared_interim(set().union(*(p.declared_elements() for p in found
                                          if _is_declared(p))))


def _builds_policy(builds, policy=None, participants=None) -> ChemistryPolicy:
    """The policy of a witness of ``builds`` (or of its ``participants``)."""
    if policy is not None:
        return policy
    if participants is None:
        _molecules, participants = b.witness_reactants(builds)
    return policy_for(_Participants(tuple(participants)))


@dataclass(frozen=True)
class _Participants:
    participants: tuple


def certificate_domain(policy: ChemistryPolicy):
    """The ``domain=`` a certificate takes under ``policy``: None (the module
    domain, certificates unchanged) for the legacy domain, the policy's
    domain otherwise."""
    return None if policy.domain is LEGACY_V1.domain else policy.domain


def _radical_ends(witness, port, providers=None) -> tuple:
    """The radical end chemistries of the port's provider."""
    return tuple(b._provider_for(witness, tuple(port), providers).RADICAL_END_CHEMISTRIES)


ROSTER_REFUSAL = "refused:roster: declared providers are not in the PS roster"


def _roster_guard(*items) -> None:
    """Refuse a declared spec, provider or all-declared witness: the named
    rosters build their own chains and never stand for a declaration."""
    for item in items:
        if item is None:
            continue
        if isinstance(item, (DeclaredSpec, dc.DeclaredChainProvider)) or \
                isinstance(item, str) and _is_declared(b.PROVIDERS.get(item)):
            raise WindowError(ROSTER_REFUSAL)
        participants = getattr(item, "participants", None)
        if participants and all(_is_declared(b.PROVIDERS.get(p.metadata.provider))
                                for p in participants):
            raise WindowError(ROSTER_REFUSAL)


def _roster_entry(make):
    def entry(provider=None):
        _roster_guard(provider)
        return make()
    return entry


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Support:
    """``atoms``: the support; ``reacting``: the centres (labels, mapped atoms,
    edits); ``sources``: (node, result, atoms) per predicate; ``dispositions``:
    (item, disposition, argument) for every non-atom-set read; ``refusal``:
    why no window can be built (None: it can)."""

    atoms: frozenset
    reacting: frozenset = frozenset()
    sources: tuple = ()
    dispositions: tuple = ()
    refusal: Optional[str] = None


@dataclass(frozen=True)
class Window:
    atoms: frozenset
    support: Support
    added: tuple = ()
    refusal: Optional[str] = None


@dataclass(frozen=True)
class Observation:
    """One observation of a witness: support, the certificate view (``viewer``:
    ``(inside, transport=None) -> hashable``), the certificate itself and
    the recertification against the base witness."""

    support: Support
    viewer: Callable = field(compare=False)
    certificate: Any = field(default=None, compare=False)
    refusal: Optional[str] = None
    recertification: Any = field(default=None, compare=False)

    def view(self, inside, transport=None):
        if transport is None:
            return self.viewer(inside)
        return self.viewer(inside, transport)


@dataclass(frozen=True)
class FamilyConstants:
    family: str
    radius: int
    reacting_distance: int
    free_label_search: bool
    bimolecular: bool = False
    directions: tuple = ()
    pair_kinds: tuple = ()    # ((direction, (kind, kind)), ...): ``pair_roster``
    classes: Optional[tuple] = None   # SiteClass, ...: ``class_roster`` (None: not computed)
    per_reactant: tuple = ()  # directions RMG rates from one tree per reactant (``per_reactant_tree``)
    joined: tuple = ()        # ((direction, unpaired delta), ...): one-product bimolecular templates
    merge: Any = None         # ClassMerge (``merge_classes``; None: every class its own representative)


@dataclass(frozen=True)
class FixedPointResult:
    status: str
    reason: str
    depths: tuple = ()
    witness: Any = field(default=None, compare=False)
    window: Optional[Window] = field(default=None, compare=False)
    observation: Optional[Observation] = field(default=None, compare=False)
    checks: tuple = field(default=(), compare=False)
    retries: int = 0
    history: tuple = field(default=(), compare=False)
    flips: tuple = field(default=(), compare=False)

    @property
    def reached(self) -> bool:
        return self.status == "fixed_point"


@dataclass(frozen=True)
class ContextResult:
    assignment: tuple
    status: str
    reason: str = ""
    witness: Any = field(default=None, compare=False)
    window: Any = field(default=None, compare=False)
    observation: Any = field(default=None, compare=False)
    fixed_point: Any = field(default=None, compare=False)


# ---------------------------------------------------------------------------
# graph helpers on witness rows
# ---------------------------------------------------------------------------
def _side_rows(witness, side):
    return witness.reactant_rows() if side == "reactants" else witness.product_rows()


def _neighbours(rows, atom):
    return [n for n, _o in rows[atom][5]]


def _pi_atoms(rows) -> set:
    return {a for a, row in rows.items()
            if row[3] > 0 or any(float(o) > 1.0 for _n, o in row[5])}


def _pi_component(rows, pi, start) -> set:
    seen, stack = {start}, [start]
    while stack:
        a = stack.pop()
        for n in _neighbours(rows, a):
            if n in pi and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def _cyclic_atoms(rows) -> set:
    """Atoms on a cycle (incident to a non-bridge bond)."""
    import networkx as nx
    graph = nx.Graph()
    graph.add_nodes_from(rows)
    for a, row in rows.items():
        for n, _o in row[5]:
            if a < n:
                graph.add_edge(a, n)
    bridges = {tuple(sorted(e)) for e in nx.bridges(graph)}
    return {a for a in rows for n in _neighbours(rows, a) if tuple(sorted((a, n))) not in bridges}


def smallest_ring(rows, atom, cyclic=None, limit: int = 64) -> Optional[tuple]:
    """Atoms of a smallest ring through ``atom`` (None if it is on none, a
    tuple longer than ``limit`` atoms is never searched: ``()`` then)."""
    cyclic = _cyclic_atoms(rows) if cyclic is None else cyclic
    if atom not in cyclic:
        return None
    best = None
    for first in _neighbours(rows, atom):
        if first not in cyclic:
            continue
        # shortest path first -> atom avoiding the bond atom-first
        parent = {first: None}
        frontier = [first]
        depth = 0
        found = False
        while frontier and depth < limit and not found:
            depth += 1
            nxt = []
            for a in frontier:
                for n in _neighbours(rows, a):
                    if a == first and n == atom:
                        continue
                    if n == atom:
                        path = [a]
                        while parent[path[-1]] is not None:
                            path.append(parent[path[-1]])
                        ring = tuple(sorted(set(path) | {atom}))
                        if best is None or len(ring) < len(best):
                            best = ring
                        found = True
                        break
                    if n not in parent:
                        parent[n] = a
                        nxt.append(n)
                if found:
                    break
            frontier = nxt
    return best if best is not None else ()


# ---------------------------------------------------------------------------
# closure (item 2)
# ---------------------------------------------------------------------------
def close(witness, support, *, ring_limit=RING_LIMIT, sides=("reactants", "products")):
    """The closed window of ``support`` (a ``Support`` or an atom set) on
    ``witness`` (see the module docstring); a refusal is returned as a
    ``Window`` with ``refusal`` set and no atoms beyond the support."""
    if not isinstance(support, Support):
        support = Support(atoms=frozenset(support))
    rows_of = {side: _side_rows(witness, side) for side in sides}
    every = set(witness.reactant_rows())
    unknown = sorted(set(support.atoms) - every)
    if unknown:
        raise WindowError(f"support atoms {unknown[:5]} are not witness atoms")
    window = set(support.atoms)
    centres = set(support.reacting) or set(support.atoms)
    info = {side: (rows, _cyclic_atoms(rows), _pi_atoms(rows)) for side, rows in rows_of.items()}
    ring_cache = {}
    added = []
    changed = True
    while changed:
        changed = False
        for side, (rows, cyclic, pi) in info.items():
            for atom in sorted(window):
                if rows[atom][0] == "H":
                    continue
                key = (side, atom)
                if key not in ring_cache:
                    ring_cache[key] = smallest_ring(rows, atom, cyclic, limit=4 * ring_limit)
                ring = ring_cache[key]
                if ring is not None:
                    if ring == () or len(ring) > ring_limit:
                        size = len(ring) if ring else f"> {4 * ring_limit}"
                        if atom in centres:
                            return Window(frozenset(support.atoms), support, tuple(added),
                                          f"refused:ring: reacting atom {atom} is in a "
                                          f"{size}-membered ring (> {ring_limit}) on the "
                                          f"{side} side")
                        return Window(frozenset(support.atoms), support, tuple(added),
                                      f"refused:large_ring: window atom {atom} is in a "
                                      f"{size}-membered ring (> {ring_limit}) on the {side} "
                                      "side (the declared domain admits phenyl rings only)")
                    new = set(ring) - window
                    if new:
                        window |= new
                        added.append(("ring", side, tuple(sorted(new))))
                        changed = True
                if atom in pi:
                    new = _pi_component(rows, pi, atom) - window
                    if new:
                        window |= new
                        added.append(("pi", side, tuple(sorted(new))))
                        changed = True
            hydrogens = {n for a in window if rows[a][0] != "H" for n in _neighbours(rows, a)
                         if rows[n][0] == "H"} - window
            if hydrogens:
                window |= hydrogens
                added.append(("H", side, tuple(sorted(hydrogens))))
                changed = True
    return Window(frozenset(window), support, tuple(added), None)


def boundary_bond_problems(witness, window) -> list:
    """Bonds leaving the window whose order is not 1 on some side (the order
    of a boundary bond must be a window property: pi systems and rings are
    closed)."""
    out = []
    atoms = window.atoms if isinstance(window, Window) else set(window)
    for side in ("reactants", "products"):
        rows = _side_rows(witness, side)
        for a in sorted(atoms):
            for n, o in rows[a][5]:
                if n not in atoms and float(o) != 1.0:
                    out.append(f"{side}: bond {a}-{n} of order {o} leaves the window")
    return out


def port_inside(witness, port_key, atoms) -> bool:
    port = witness.port(port_key)
    return port.atom_id in atoms or any(c in atoms for c in port.cap_atom_ids)


# ---------------------------------------------------------------------------
# window marked graphs
# ---------------------------------------------------------------------------
def _port_atoms(witness) -> dict:
    out = {}
    for participant in witness.participants:
        for port in participant.metadata.ports:
            out[port.atom_id] = (participant.index, port.port_id)
            for cap in port.cap_atom_ids:
                out[cap] = (participant.index, port.port_id)
    return out


def _roles(witness) -> dict:
    return {a: p.index for p in witness.participants for a in p.atom_ids}


def window_labels(witness, atoms) -> tuple:
    """(node labels {atom: label}, edge labels {(a, b): label}) of the window's
    marked graph: attributes on both sides, participant role, reacting mark
    (the edit support), physical-end mark, and per exterior bond its orders
    and the ports its far side reaches (orientation)."""
    atoms = set(atoms)
    before, after = witness.reactant_rows(), witness.product_rows()
    support = b.edit_support(witness)
    roles = _roles(witness)
    ports = _port_atoms(witness)
    physical = {}
    for participant in witness.participants:
        for port in participant.metadata.ports:
            if not port.artificial:
                physical[port.atom_id] = (participant.index, port.port_id, port.tag)
    bonds_before, bonds_after = witness.reactant_bonds(), witness.product_bonds()
    adjacency = defaultdict(set)
    for pair in set(bonds_before) | set(bonds_after):
        adjacency[pair[0]].add(pair[1])
        adjacency[pair[1]].add(pair[0])
    reach_cache = {}

    def reach(start):
        # port keys reachable from ``start`` outside the window (either side)
        if start in reach_cache:
            return reach_cache[start]
        seen, stack, found = {start}, [start], set()
        while stack:
            a = stack.pop()
            if a in ports:
                found.add(ports[a])
            for n in adjacency[a]:
                if n not in atoms and n not in seen:
                    seen.add(n)
                    stack.append(n)
        value = tuple(sorted(found))
        for a in seen:
            reach_cache[a] = value
        return value

    nodes = {}
    for a in sorted(atoms):
        exterior = []
        for n in sorted(adjacency[a]):
            if n in atoms:
                continue
            pair = (min(a, n), max(a, n))
            exterior.append((bonds_before.get(pair, 0.0), bonds_after.get(pair, 0.0), reach(n)))
        nodes[a] = (b._attrs(before[a]), b._attrs(after[a]), roles.get(a, -1),
                    a in support, physical.get(a), tuple(sorted(exterior, key=repr)))
    edges = {}
    for a in atoms:
        for n in adjacency[a]:
            if n in atoms and a < n:
                edges[(a, n)] = (bonds_before.get((a, n), 0.0), bonds_after.get((a, n), 0.0))
    return nodes, edges


def same_window(first, atoms_first, second, atoms_second) -> bool:
    """Identity correspondence (surviving atoms keep their ids)."""
    if set(atoms_first) != set(atoms_second):
        return False
    return window_labels(first, atoms_first) == window_labels(second, atoms_second)


def _collapsed_graph(witness, atoms):
    """networkx graph of the window: heavy atoms and reacting H as nodes,
    the other H folded into their heavy neighbour's label; reacting atoms
    carry their id (marks fix them)."""
    import networkx as nx
    nodes, edges = window_labels(witness, atoms)
    rows = witness.reactant_rows()
    support = b.edit_support(witness)
    kept = {a for a in nodes if rows[a][0] != "H" or a in support}
    folded = defaultdict(list)
    for a in nodes:
        if a in kept:
            continue
        heavy = [n for n, _o in rows[a][5] if n in kept]
        folded[heavy[0] if heavy else None].append(repr(nodes[a]))
    graph = nx.Graph()
    for a in kept:
        mark = f"r{a}" if a in support else ""
        graph.add_node(a, label=repr((nodes[a], mark, tuple(sorted(folded.get(a, ()))))))
    graph.graph["orphans"] = tuple(sorted(folded.get(None, ())))
    for (x, y), label in edges.items():
        if x in kept and y in kept:
            graph.add_edge(x, y, label=repr(label))
    return graph


def window_isomorphisms(first, atoms_first, second, atoms_second, limit=ISOMORPHISM_LIMIT):
    """Mark-preserving isomorphisms (second -> first ids over heavy atoms and
    reacting H) of two windows' marked graphs (at most ``limit``)."""
    from networkx.algorithms import isomorphism
    g1, g2 = _collapsed_graph(first, atoms_first), _collapsed_graph(second, atoms_second)
    if g1.graph["orphans"] != g2.graph["orphans"] or len(g1) != len(g2):
        return []
    matcher = isomorphism.GraphMatcher(
        g2, g1, node_match=lambda x, y: x["label"] == y["label"],
        edge_match=lambda x, y: x["label"] == y["label"])
    out = []
    for mapping in matcher.isomorphisms_iter():
        out.append(dict(mapping))
        if len(out) >= limit:
            break
    return out


# ---------------------------------------------------------------------------
# certificate support and view (item 1)
# ---------------------------------------------------------------------------
def _label_set(mapping) -> frozenset:
    return frozenset((lab, a) for _s, lab, a in mapping["labels"])


def member_entries(cert, wanted) -> list:
    """(channel, member) of every certified member one of whose mappings has
    a label set in ``wanted``."""
    out = []
    for ch in cert.get("channels", []):
        for mem in ch.get("members", []):
            if any(_label_set(mp) in wanted for mp in mem["mappings"]):
                out.append((ch, mem))
    return out


def mapping_atoms(mapping) -> set:
    atoms = {a for _s, _l, a in mapping["labels"]}
    atoms |= {a for _s, _i, _j, a in mapping.get("mapped", [])}
    edits = mapping.get("edits") or {}
    atoms |= {a for a, _b, _c in edits.get("atoms", [])}
    for a, c, _o1, _o2 in edits.get("bonds", []):
        atoms |= {a, c}
    return atoms


class SearchClassifier:
    """Whether an extension can add a candidate to a free-label search:
    tests the RMG group atoms of the label against every atom of an interior
    repeat unit of each provider (RMG atom types)."""

    def __init__(self, family, providers: Optional[Sequence] = None):
        self.family = family
        self.providers = tuple(provider_object(p) for p in (
            b.LEGACY_PROVIDERS if providers is None else providers))
        self.fresh = []
        for provider in self.providers:
            self.fresh.extend(_interior_unit_atoms(provider))
        self._cache = {}

    def extensible(self, node, label) -> bool:
        key = (node, label)
        if key not in self._cache:
            from rmgpy.kmc import channels as chm
            entry = self.family.groups.entries.get(node)
            if entry is None:
                self._cache[key] = True    # unknown semantics: never claimed bounded
            else:
                hit = False
                for group in chm._group_atoms(self.family.groups, entry):
                    centres = group.get_all_labeled_atoms().get(label)
                    if centres is None:
                        continue
                    gatoms = centres if isinstance(centres, list) else [centres]
                    if any(atom.is_specific_case_of(g) for atom in self.fresh for g in gatoms):
                        hit = True
                        break
                self._cache[key] = hit
        return self._cache[key]


_INTERIOR = {}
_FACTS = {}
FACTS_UNDECIDED = "undecided"


def _interior_unit_atoms(provider) -> list:
    """RMG atoms (atom types perceived) of an interior repeat of a provider
    (an object or a registered name) and their H (``interior_atoms_rmg``),
    cached by provider name."""
    provider = provider_object(provider)
    if provider.name not in _INTERIOR:
        _INTERIOR[provider.name] = list(provider.interior_atoms_rmg())
    return _INTERIOR[provider.name]


def _fragment_facts(providers=None) -> dict:
    """Per provider (objects; default the named pair): whether every extension
    fragment of its probes (``probe_ports``) is closed shell and the net
    element delta of one step (fresh atoms minus the deleted caps).  A probe
    that raises leaves that provider's facts undecided
    (``{FACTS_UNDECIDED: reason}``).  Cached by the tuple of provider
    names."""
    providers = tuple(provider_object(p) for p in (
        b.LEGACY_PROVIDERS if providers is None else providers))
    key = tuple(p.name for p in providers)
    if key not in _FACTS:
        out = {}
        for provider in providers:
            facts = {}
            try:
                for port, coordinate in provider.probe_ports():
                    fragment = provider.fragment(port, coordinate)
                    delta = defaultdict(int)
                    for atom in fragment.atoms:
                        delta[atom[0]] += 1
                    delta["H"] -= len(port.cap_atom_ids)
                    facts[port.tag] = {"closed_shell": all(atom[3] == 0 for atom in fragment.atoms),
                                       "delta": dict(delta)}
            except Exception as error:   # fail closed: the facts are not decided
                facts = {FACTS_UNDECIDED: f"{type(error).__name__}: {error}"}
            out[provider.name] = facts
        _FACTS[key] = out
    return _FACTS[key]


def certificate_support(cert, wanted, *, classifier: Optional[SearchClassifier] = None,
                        fragment_facts=None, include_failed: bool = True,
                        growable=None) -> Support:
    """The support of the member of ``cert`` selected by ``wanted`` (label
    sets), with the disposition of every non-atom-set read (see the module
    docstring).  ``growable``: the witness atoms of participants an
    extension can grow (an artificial port; ``growable_atoms``); an
    ``element_count`` guard read on a certificate molecule without such an
    atom is extension invariant (no extension adds to that molecule).  None
    (the default): every molecule takes every provider's step delta."""
    entries = member_entries(cert, wanted)
    if not entries:
        return Support(frozenset(), refusal="refused:member: no certified member carries the "
                                            "wanted labels")
    if len({id(ch) for ch, _m in entries}) > 1:
        return Support(frozenset(), refusal=f"refused:member: the member is in "
                                            f"{len(entries)} channels")
    facts = fragment_facts if fragment_facts is not None else _fragment_facts()
    undecided = sorted(name for name, p in facts.items() if FACTS_UNDECIDED in p)
    decided = [p for p in facts.values() if FACTS_UNDECIDED not in p]
    closed_shell = all(f["closed_shell"] for p in decided for f in p.values())
    deltas = defaultdict(set)
    for p in decided:
        for f in p.values():
            for el, d in f["delta"].items():
                deltas[el].add(d)
    molecule_of, still = {}, set()
    if growable is not None:
        for k, ids in enumerate(cert["witness"]["molecules"]):
            molecule_of.update((int(a), k) for a in ids)
            if not {int(a) for a in ids} & set(growable):
                still.add(k)
    reacting, atoms, sources = set(), set(), []
    dispositions = {}
    refusal = None
    for ch, mem in entries:
        for mp in mem["mappings"]:
            reacting |= mapping_atoms(mp)
    for ch, mem in entries:
        for mp in mem["mappings"]:
            if mp.get("descent") is None:
                refusal = refusal or f"refused:incomplete: mapping without descent: {mp.get('gap')}"
                continue
            for d in cert["descents"][mp["descent"]]:
                if d.get("skipped"):
                    continue
                for p in _iter_predicates(d):
                    if p.get("logic"):
                        continue
                    if p["result"] == "failed" and not include_failed:
                        continue
                    pa = set((p.get("anchors") or {}).values())
                    pa |= {i for i, _d in p.get("ball", [])}
                    emb = p.get("embeddings") or {}
                    pa |= set(emb.get("images", []))
                    pa |= set(p.get("ignored", []))
                    for lab in p.get("free_labels", []):
                        cands = set((p.get("candidates") or {}).get(lab, []))
                        key = ("search", p["node"], lab)
                        if classifier is None or classifier.extensible(p["node"], lab):
                            dispositions[key] = ("refused", "an extension can add candidates "
                                                 "of this free label (unbounded search)")
                            refusal = refusal or (f"refused:unbounded_search: {p['node']} "
                                                  f"label {lab}")
                        elif not cands <= reacting:
                            dispositions[key] = ("whole_species", "a candidate outside the "
                                                 f"reacting atoms: {sorted(cands - reacting)}")
                            refusal = refusal or (f"refused:whole_species_search: {p['node']} "
                                                  f"label {lab} candidates "
                                                  f"{sorted(cands - reacting)[:4]}")
                        else:
                            dispositions[key] = ("bounded_by_domain", "no interior repeat-unit "
                                                 "atom matches the free group atom (no extension "
                                                 "adds a candidate); every candidate is a "
                                                 "reacting atom")
                    for name, g in sorted((p.get("guards") or {}).items()):
                        if not g.get("read"):
                            continue
                        if name in ("multiplicity", "radical_count"):
                            if undecided:
                                dispositions[("guard", name)] = (
                                    "refused", f"fragment facts undecided for {undecided[0]}")
                                refusal = refusal or \
                                    f"refused:fragment_facts_undecided: {undecided[0]}"
                            elif not closed_shell:
                                dispositions[("guard", name)] = ("refused", "a provider "
                                                                 "fragment is not closed shell")
                                refusal = refusal or f"refused:guard: {name}"
                            else:
                                dispositions[("guard", name)] = (
                                    "context_condition", "whole-species unpaired-electron count "
                                    "(declared maximum multiplicity); extension invariant "
                                    "(closed-shell fragments); varied by the contexts' radical "
                                    "inventory; its value is in the certificate view")
                        elif name == "element_count":
                            minimum = g.get("group_minimum") or {}
                            have = g.get("structure") or {}
                            read = {molecule_of[a] for a in pa if a in molecule_of}
                            invariant = bool(read) and read <= still
                            for el, need in sorted(minimum.items()):
                                count = have.get(el, 0)
                                ds = {0} if invariant else deltas.get(el, {0})
                                key = ("guard", f"element_count:{el}")
                                if undecided:
                                    dispositions[key] = ("refused", "fragment facts undecided "
                                                         f"for {undecided[0]}")
                                    refusal = refusal or \
                                        f"refused:fragment_facts_undecided: {undecided[0]}"
                                elif any(x < 0 for x in ds):
                                    dispositions[key] = ("refused", "an extension removes "
                                                         f"{el} atoms")
                                    refusal = refusal or f"refused:guard: element_count {el}"
                                elif count >= need:
                                    dispositions[key] = ("monotone_satisfied", "satisfied; "
                                                         "extension never lowers the count")
                                elif ds == {0}:
                                    dispositions.setdefault(key, (
                                        "invariant_failed", "failed; extension never adds "
                                        f"{el}"))
                                else:
                                    dispositions[key] = ("refused", f"failed ({count} < {need}); "
                                                         "an extension could satisfy it")
                                    refusal = refusal or f"refused:guard: element_count {el}"
                    sources.append((p["node"], p["result"], frozenset(pa)))
                    atoms |= pa
    atoms |= reacting
    for ch, _mem in entries:
        dispositions[("global", "channel_kind")] = (
            "whole_species_identity", f"{ch.get('kind')}: kept vs synthetic:unchanged_species "
            "(products equal to the reactant species) and the identical-reactant pairing are "
            "decided by whole-species identity; both kinds are rated by the family's own "
            "estimate for the template (channels.channel_kinetics), which the view carries")
    for guard in cert.get("global_guards", []):
        name = guard["guard"]
        if name == "multiplicity_convention":
            dispositions[("global", name)] = ("declared", "domain.MULTIPLICITY_CONVENTION, "
                                              "checked on every species by domain.check")
        elif name == "inRing":
            dispositions[("global", name)] = ("closure", "every ring <= RING_LIMIT through a "
                                              "window atom is in the window")
        elif name == "depository_whole_species":
            for ch, _mem in entries:
                if ch.get("whole_species"):
                    dispositions[("global", name)] = ("whole_species", "a depository entry "
                                                      "matched the whole species")
                    refusal = refusal or "refused:whole_species: depository match"
                else:
                    verdict = (ch.get("depository_exclusion") or {}).get("verdict", "")
                    if verdict.startswith("refused") or not verdict:
                        dispositions[("global", name)] = ("refused", f"depository exclusion "
                                                          f"{verdict or 'missing'}")
                        refusal = refusal or f"refused:depository_exclusion: {verdict}"
                    else:
                        dispositions[("global", name)] = ("bounded_by_domain",
                                                          f"depository exclusion: {verdict}")
    return Support(frozenset(atoms), frozenset(reacting), tuple(sources),
                   tuple(sorted((k, v[0], v[1]) for k, v in dispositions.items())), refusal)


def _iter_predicates(descent):
    from rmgpy.kmc import certificates as cm
    return cm.iter_descent_predicates(descent)


def growable_atoms(witness) -> frozenset:
    """The atoms of the participants of ``witness`` with an artificial port
    (the molecules an extension can grow)."""
    ports = {key[0] for key in b.vector_ports(witness)}
    return frozenset(a for p in witness.participants if p.index in ports for a in p.atom_ids)


def certificate_refusal(cert) -> Optional[str]:
    """Why a certificate cannot carry a window (None: it can)."""
    if cert is None:
        return "refused:certificate: none"
    if not cert.get("complete"):
        return f"refused:incomplete: {(cert.get('gaps') or ['?'])[0]}"[:300]
    if not cert.get("eligible"):
        why = cert.get("input_refusals") or [
            (ch.get("eligibility") or {}).get("disposition") for ch in cert["channels"]
            if (ch.get("eligibility") or {}).get("disposition") != "eligible"]
        return f"refused:domain: {why[:2]}"[:300]
    return None


def member_refusal(cert, wanted) -> Optional[str]:
    for ch, _mem in member_entries(cert, wanted):
        if ch.get("support") != "supported":
            return f"refused:{ch.get('support')}"
        if not ch.get("publishable", True):
            return f"refused:unpublishable: {(ch.get('eligibility') or {}).get('disposition')}"
    return None


def certificate_view(cert, wanted, inside, transport=None, *, transport_reacting=False):
    """The member's certificate, correspondence-aware, restricted to the
    window ``inside``: per mapping (identity with ids through the map) its
    channel selection and status, descent signature (every predicate tried,
    verdicts, choices, leaves) and every predicate's local support (anchors,
    candidates, ball, embeddings, ignored atoms) and guard verdicts, and the
    atom-table rows of window atoms.  Non-reacting H are named by their heavy
    neighbour; atoms outside the window are ``out``.  Reacting atoms keep
    their ids (contexts of one witness) unless ``transport_reacting`` (class
    equivalence across witnesses: they go through ``transport`` too)."""
    from rmgpy.kmc import certificates as cm
    inside = set(inside)
    entries = member_entries(cert, wanted)
    reacting = set()
    for _ch, mem in entries:
        for mp in mem["mappings"]:
            reacting |= mapping_atoms(mp)
    witness_atoms = {row[0]: row for row in cert["witness"]["atoms"]}
    base = (lambda a: a) if transport is None else (lambda a: transport.get(a, ("unmapped", a)))

    def t(a):
        if a in reacting:
            return ("r", base(a) if transport_reacting else a)
        if a not in inside:
            return "out"
        row = witness_atoms.get(a)
        if row is not None and row[1] == "H" and row[3]:
            return ("H", t(row[3][0]))
        return base(a)

    def tl(values):
        return tuple(sorted((t(v) for v in values), key=repr))

    out = []
    for ch, mem in entries:
        # the channel kind (kept | synthetic:unchanged_species) is a
        # whole-species identity verdict, not a read of the selection: it is
        # a disposition (``certificate_support``), not part of the view
        selection = list(cm._channel_selection(ch))
        selection[1] = None
        sel = tuple(selection) + (
            ch.get("support"), (ch.get("eligibility") or {}).get("disposition"),
            ch.get("whole_species"), (ch.get("depository_exclusion") or {}).get("verdict"))
        for mp in mem["mappings"]:
            descents = cert["descents"].get(mp.get("descent"), []) if mp.get("descent") else []
            preds = []
            rows = []
            for d in descents:
                if d.get("skipped"):
                    continue
                for p in _iter_predicates(d):
                    if p.get("logic"):
                        preds.append((p["node"], p["result"], "logic"))
                        continue
                    emb = p.get("embeddings") or {}
                    guards = []
                    for name, g in sorted((p.get("guards") or {}).items()):
                        if name == "element_count":
                            have = g.get("structure") or {}
                            guards.append((name, g.get("read"), tuple(sorted(
                                (el, have.get(el, 0) >= need)
                                for el, need in (g.get("group_minimum") or {}).items()))))
                        else:
                            guards.append((name, g.get("read"), g.get("structure")))
                    preds.append((p["node"], p["result"], p.get("anchoring"),
                                  tuple(sorted((lab, t(a)) for lab, a in
                                               (p.get("anchors") or {}).items())),
                                  tuple(sorted((lab, tl(ids)) for lab, ids in
                                               (p.get("candidates") or {}).items())),
                                  tuple(sorted(((t(i), dd) for i, dd in p.get("ball", [])),
                                               key=repr)),
                                  (emb.get("count"),
                                   tuple(sorted((int(g), tl(ids))
                                                for g, ids in emb.get("by_group", []))),
                                   tl(emb.get("images", []))),
                                  tuple(guards), tl(p.get("ignored", []))))
                for row in d.get("atoms", []):
                    if row[0] in inside:
                        rows.append((t(row[0]),) + tuple(row[1:8]) + (
                            tuple(sorted(((t(n), o) for n, o in row[8]), key=repr)),))
            key = cm.mapping_key(mp, t)
            out.append((key, sel, cm.descent_signature(descents), tuple(preds),
                        tuple(sorted(rows, key=repr))))
    return tuple(sorted(out, key=repr))


# ---------------------------------------------------------------------------
# observers
# ---------------------------------------------------------------------------
def connectivity_change(family, labels) -> frozenset:
    """The bonds the FORWARD (family-template) recipe makes or breaks, as atom
    pairs, for a template-side label set ``{(label, atom id), ...}`` (a label
    may name two atoms: R_Recombination's two ``*``).  Two mappings with the
    same mapped outcome make and break the same bonds (the outcome fixes the
    product connectivity; resonance moves only bond orders and electrons)."""
    atoms = defaultdict(list)
    for label, atom in labels:
        atoms[label].append(atom)
    out = set()
    for action in family.forward_recipe.actions:
        if action[0] not in ("FORM_BOND", "BREAK_BOND"):
            continue
        first, second = action[1], action[3]
        if first == second:
            pairs = [tuple(sorted(atoms[first]))] if len(atoms[first]) == 2 else []
        else:
            pairs = [tuple(sorted((x, y))) for x in atoms[first] for y in atoms[second]]
        out.update((action[0], pair) for pair in pairs)
    return frozenset(out)


def _candidate_labels(candidate) -> frozenset:
    return frozenset((lab, a) for _s, lab, a in candidate.labels)


class ScopedCertifier:
    """A ``certificates.Certifier`` whose forward channel sets keep only the
    mappings that make and break the member's bonds (``connectivity_change``):
    every mapping of the member (all slot orders, all resonance variants) and
    every mapping of any other pathway with the same mapped outcome (pathway
    ambiguity), none of the unrelated ones.  The member's mappings, descents,
    selection and channel status are those of the full certificate (checked:
    ``windowsTest``); the channel's other members (other bonds) are other
    windows.  Shares every other cache with the full certifier."""

    def __new__(cls, base, family, change):
        from rmgpy.kmc import certificates as cm

        class _Scoped(cm.Certifier):
            def channel_set(self, family_, reactants):
                from rmgpy.kmc import channels as chm
                from rmgpy.kmc import rmg_moves as rm
                if family_.label != self._scope_family:
                    return cm.Certifier.channel_set(self, family_, reactants)
                roots = [cm.aromatic_root(m) for m in reactants]
                key = ("scoped", family_.label, tuple(sorted(self._scope_change)),
                       tuple(rm._signature(rm.molecule_rows(m)) for m in roots))
                cs = self._channel_sets.get(key)
                if cs is None:
                    cs = chm.build_channels(
                        family_, roots, "forward", resonance=self.resonance,
                        select=lambda c: connectivity_change(
                            family_, _candidate_labels(c)) == self._scope_change)
                    self._channel_sets[key] = cs
                return cs

        scoped = _Scoped.__new__(_Scoped)
        scoped.__dict__.update(base.__dict__)
        scoped._scope_family = family.label
        scoped._scope_change = frozenset(change)
        return scoped


def scoped_certifier(certifier, family, wanted):
    """A ``ScopedCertifier`` for the member whose template-side label sets are
    ``wanted`` (all must share one connectivity change)."""
    changes = {connectivity_change(family, labels) for labels in wanted}
    if len(changes) != 1:
        raise WindowError(f"the member's mappings make different bonds: {changes}")
    return ScopedCertifier(certifier, family, changes.pop())


def template_molecules(witness, direction: str):
    """RMG molecules of the family-template side of ``witness``: its reactants
    for a forward record, its products for a reverse record (components by
    least atom id; atom ids kept)."""
    from rmgpy.kmc import rmg_moves as rm
    rows = witness.reactant_rows() if direction == "forward" else witness.product_rows()
    return [rm.molecule_from_rows(part) for part in rm.split_rows(dict(rows))]


class MemberObserver:
    """Observe a member on a witness: certify the family-template side (the
    first observation is the base; later ones are ``recertify_certificate``
    against it through the identity map of surviving atoms), then support,
    dispositions and the certificate view."""

    def __init__(self, certifier, family, direction, wanted, *, base_certificate=None,
                 classifier=None, include_failed=True, scoped=True, providers=None,
                 policy=None):
        self.wanted = frozenset(wanted)
        self.certifier = scoped_certifier(certifier, family, self.wanted) if scoped else \
            certifier
        self.family = family
        self.direction = direction
        self.base = base_certificate
        self.base_atoms = None if base_certificate is None else \
            {row[0] for row in base_certificate["witness"]["atoms"]}
        # resolved from the first observed witness's participants
        self.classifier = classifier
        self.facts = None
        self.repeat = None
        self.providers = None
        self.registry = providers
        self.include_failed = include_failed
        self.calls = 0
        self.flips = []
        # resolved from the first observed witness unless given
        self.policy = policy
        self.domain = None

    def __call__(self, witness) -> Observation:
        self.calls += 1
        if self.providers is None:
            self.providers = observation_providers(witness, self.registry)
            if self.classifier is None:
                self.classifier = SearchClassifier(self.family, providers=self.providers)
            self.facts = _fragment_facts(self.providers)
            self.repeat = growth_unit(witness, self.registry)
            if self.policy is None:
                self.policy = policy_for(witness, self.registry)
            self.domain = certificate_domain(self.policy)
        extra = {} if self.repeat is None else {"repeat": self.repeat}
        if self.domain is not None:
            extra["domain"] = self.domain
        mols = template_molecules(witness, self.direction)
        recert = None
        if self.base is None:
            cert = self.certifier.certify(self.family, mols, member_labels=self.wanted, **extra)
            self.base = cert
            self.base_atoms = {row[0] for row in cert["witness"]["atoms"]}
        else:
            present = {a.id for m in mols for a in m.atoms}
            amap = {a: a for a in self.base_atoms if a in present}
            recert = self.certifier.recertify_certificate(self.family, self.base, mols, amap,
                                                          **extra)
            cert = recert["certificate"]
            for f in recert.get("flips", []):
                self.flips.append((dict(witness.extension_counts), f["node"], f["old"], f["new"]))
            if recert.get("refusals"):
                return Observation(Support(frozenset()), lambda inside, t=None: None, cert,
                                   "refused:recertification: " + recert["refusals"][0][:240],
                                   recert)
        refusal = certificate_refusal(cert) or member_refusal(cert, self.wanted)
        support = certificate_support(cert, self.wanted, classifier=self.classifier,
                                      fragment_facts=self.facts,
                                      include_failed=self.include_failed,
                                      growable=growable_atoms(witness))
        refusal = refusal or support.refusal
        missing = sorted(set(support.atoms) - set(witness.reactant_rows()))
        if missing and not refusal:
            refusal = f"refused:transport: support atoms {missing[:4]} are not witness atoms"
        wanted = self.wanted

        def viewer(inside, transport=None, cert=cert):
            return certificate_view(cert, wanted, inside, transport)
        return Observation(support, viewer, cert, refusal, recert)


# ---------------------------------------------------------------------------
# fixed point (item 3)
# ---------------------------------------------------------------------------
def _key(witness):
    return tuple(sorted((k, v) for k, v in witness.extension_counts if v))


def check_groups(ports) -> tuple:
    """The checks of one round: (ports, step, joint) for every artificial port
    alone at every step of ``CHECK_STEPS``, then all ports jointly at every
    step.  No option removes any of them."""
    ports = tuple(ports)
    out = [((p,), step, False) for p in ports for step in CHECK_STEPS]
    if ports:
        out += [(ports, step, True) for step in CHECK_STEPS]
    return tuple(out)


def fixed_point(witness, observe, providers=None, *, ceiling=None, port_limit=PORT_DEPTH_LIMIT,
                retries=FIXED_POINT_RETRIES, closure=None) -> FixedPointResult:
    """The recertified fixed point of ``witness`` under ``observe`` (see the
    module docstring): every round checks every artificial port alone at +1
    and +2 repeats and all ports jointly at +1 and +2 (``check_groups``)."""
    closure = closure or close
    limit = b._ceiling(ceiling)
    ports = b.vector_ports(witness)
    cache = {}
    history = []

    def look(w):
        key = _key(w)
        if key not in cache:
            obs = observe(w)
            win = None if obs.refusal else closure(w, obs.support)
            cache[key] = (obs, win)
            history.append((dict(key), None if win is None else len(win.atoms),
                            obs.refusal or (win.refusal if win else None)))
        return cache[key]

    def depths(w):
        return tuple((p, w.depth(p)) for p in ports)

    def stop(status, reason, w=None, win=None, obs=None, checks=(), used=0):
        return FixedPointResult(status, reason, depths(w or witness), w, win, obs,
                                tuple(checks), used, tuple(history),
                                tuple(getattr(observe, "flips", ())))

    def grow(w, keys, steps=1, slack=0):
        # ``slack``: a check witness may reach one repeat past the per-port
        # limit (its +2 step), the fixed point itself never
        for p in keys:
            if w.depth(p) + steps > port_limit + slack:
                return None, stop("no_fixed_point", f"port {p} reached the per-port limit "
                                  f"{port_limit} without a fixed point", w)
        try:
            for _step in range(steps):
                for p in keys:
                    w, _record = b.extend(w, p, providers, ceiling=limit)
        except b.CeilingExceeded as error:
            return None, stop("ceiling_exceeded", f"{error} (never truncated)", w)
        except b.ExtensionRefused as error:
            if "reacting atom" in str(error):
                return None, stop("refused:reacting_port", f"a reacting atom sits on an "
                                  f"artificial port (a cut-point event: its long-chain and "
                                  f"physical-end analogues are other witnesses): {error}", w)
            return None, stop("refused:extension", f"an artificial port in the window cannot "
                              f"be extended: {error}", w)
        return w, None

    current = witness
    obs, win = look(current)
    if obs.refusal or win.refusal:
        return stop(obs.refusal or win.refusal, "refused", current, win, obs)
    rounds = []
    used = 0
    while True:
        inside = [p for p in ports if port_inside(current, p, win.atoms)]
        while inside:
            current, failure = grow(current, inside)
            if failure:
                return failure
            obs, win = look(current)
            if obs.refusal or win.refusal:
                return stop(obs.refusal or win.refusal, "refused on extension", current, win,
                            obs, rounds, used)
            inside = [p for p in ports if port_inside(current, p, win.atoms)]
        view = obs.view(win.atoms)
        checks = []
        for group, step, joint in check_groups(ports):
            other, failure = grow(current, group, step, slack=step - 1)
            if failure:
                if failure.status == "ceiling_exceeded":
                    return failure
                return stop(failure.status, failure.reason, current, None, obs, rounds, used)
            o2, w2 = look(other)
            if o2.refusal or w2.refusal:
                return stop("certificate_unstable", f"+{step} at {list(group)}: "
                            f"{o2.refusal or w2.refusal}", current, None, obs, rounds, used)
            checks.append({"ports": tuple(group), "step": step, "joint": joint,
                           "same_window": same_window(current, win.atoms, other, w2.atoms),
                           "same_certificate": view == o2.view(w2.atoms),
                           "window_atoms": len(w2.atoms)})
        rounds.append(tuple(checks))
        failed = [c for c in checks if not (c["same_window"] and c["same_certificate"])]
        if not failed:
            return stop("fixed_point", "", current, win, obs, rounds, used)
        if used >= retries:
            unstable = any(not c["same_certificate"] for c in failed)
            return stop("certificate_unstable" if unstable else "no_fixed_point",
                        f"checks still fail after {used} re-extensions: "
                        f"{[(list(c['ports']), c['step']) for c in failed]}", current, None,
                        obs, rounds, used)
        used += 1
        # failed independent checks extend their own port; a failure seen only
        # jointly extends every port
        alone = sorted({p for c in failed if not c["joint"] for p in c["ports"]})
        current, failure = grow(current, alone or sorted(ports))
        if failure:
            return failure
        obs, win = look(current)
        if obs.refusal or win.refusal:
            return stop(obs.refusal or win.refusal, "refused on re-extension", current, win, obs,
                        rounds, used)


# ---------------------------------------------------------------------------
# contexts (item 5)
# ---------------------------------------------------------------------------
def context_assignments(witness, near_ports=None, providers=None) -> tuple:
    """Assignments ((port, end chemistry or None), ...) over the artificial
    ports (vector order).  The full product over ALL ports, across
    participants, whatever the certificate reads: a near port (``near_ports``;
    default every artificial port) takes the artificial cut or any declared
    end; a far port the artificial cut or a radical end (the radical
    inventory is read whole-species).  Then, per far port, one single per
    non-radical end: the check of the far reduction (``far_port_verdicts``).
    The first assignment is all artificial."""
    ports = list(b.vector_ports(witness))
    chem = {p: tuple(b.end_chemistries(witness, p, providers)) for p in ports}
    radical = {p: _radical_ends(witness, p, providers) for p in ports}
    near = set(ports) if near_ports is None else set(near_ports) & set(ports)
    choices = [(None,) + chem[p] if p in near else
               (None,) + tuple(e for e in chem[p] if e in radical[p]) for p in ports]
    seen, out = set(), []

    def add(assignment):
        if assignment not in seen:
            seen.add(assignment)
            out.append(assignment)
    for combo in product(*choices):
        add(tuple(zip(ports, combo)))
    for p in ports:
        if p in near:
            continue
        for e in chem[p]:
            if e not in radical[p]:
                add(tuple((q, e if q == p else None) for q in ports))
    # the joint check of the far reduction: every combination of non-radical
    # ends over all far ports of one participant, and over all far ports
    for group in far_port_groups(ports, near):
        options = [[e for e in chem[p] if e not in radical[p]] for p in group]
        for combo in product(*options):
            choice = dict(zip(group, combo))
            add(tuple((q, choice.get(q)) for q in ports))
    return tuple(out)


def far_port_groups(ports, near) -> tuple:
    """The far-port groups whose non-radical ends are checked jointly: the far
    ports of each participant (two or more) and all far ports (two or more,
    across participants)."""
    far = [p for p in ports if p not in set(near)]
    by_part = defaultdict(list)
    for p in far:
        by_part[p[0]].append(p)
    groups = [tuple(g) for _k, g in sorted(by_part.items()) if len(g) > 1]
    if len(far) > 1 and tuple(far) not in groups:
        groups.append(tuple(far))
    return tuple(groups)


def far_port_verdicts(results, far, providers=None, *, witness=None) -> dict:
    """Whether each far non-radical context merges with the all-artificial
    context ``results[0]`` (``contexts_merge``: equal closed windows as marked
    graphs and equal certificate views): singles keyed (far port, end), joint
    contexts (two or more far ports set, all non-radical: single-end merging
    does not imply that combinations merge) keyed ((port, end), ...).  A
    refused or unmerged context is False: the far reduction does not hold
    there.  The radical ends are the port provider's
    (``RADICAL_END_CHEMISTRIES``), resolved on ``witness`` (default: the
    first context witness; every context has the same participants)."""
    far = set(far)
    out = {}
    if witness is None:     # the port providers: any context's witness (same participants)
        witness = next((r.witness for r in results if r.witness is not None), None)
    radical = {}
    for r in results[1:]:
        set_ = tuple((p, e) for p, e in r.assignment if e is not None)
        for p, _e in set_:
            if p not in radical:
                if witness is None:
                    raise WindowError("refused:far_ports: no context witness resolves the "
                                      "port providers")
                radical[p] = _radical_ends(witness, p, providers)
        if not set_ or any(p not in far or e in radical[p] for p, e in set_):
            continue
        out[set_[0] if len(set_) == 1 else set_] = contexts_merge(results[0], r)
    return out


def _verdict_ports(key) -> tuple:
    return (key[0],) if isinstance(key[1], str) else tuple(p for p, _e in key)


def far_port_promotions(verdicts) -> tuple:
    """Far ports whose reduction failed its check (alone or jointly): they
    become near."""
    return tuple(sorted({p for key, ok in verdicts.items() if not ok
                         for p in _verdict_ports(key)}))


def reconstruct(witness, assignment, providers=None):
    """The witness with every assigned end built (``boundaries.reconstruct_end``)."""
    for port, end in assignment:
        if end is not None:
            witness, _record = b.reconstruct_end(witness, port, end, providers)
    return witness


def contexts_merge(first: ContextResult, second: ContextResult) -> bool:
    if first.status != "fixed_point" or second.status != "fixed_point":
        return False
    w1, w2 = first.window, second.window
    if set(w1.atoms) == set(w2.atoms) and \
            window_labels(first.witness, w1.atoms) == window_labels(second.witness, w2.atoms):
        if first.observation.view(w1.atoms) == second.observation.view(w2.atoms):
            return True
    for mapping in window_isomorphisms(first.witness, w1.atoms, second.witness, w2.atoms):
        if first.observation.view(w1.atoms, {a: a for a in w1.atoms}) == \
                second.observation.view(w2.atoms, mapping):
            return True
    return False


def merge_contexts(results) -> tuple:
    """(representative index, merged index) pairs: a context merges into the
    first earlier class representative whose closed window is equal as a
    marked graph and whose certificate view is equal under that
    correspondence."""
    representatives, merges = [], []
    for j, result in enumerate(results):
        if result.status != "fixed_point":
            continue
        for i in representatives:
            if contexts_merge(results[i], result):
                merges.append((i, j))
                break
        else:
            representatives.append(j)
    return tuple(merges)


# ---------------------------------------------------------------------------
# span bound (item 4) and ladder (item 6)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class SpanBound:
    family: str
    radius: int
    reacting_distance: int
    reach: int
    backbone_atoms: int
    repeats: int
    L: int
    formula: str
    argument: str


SPAN_ARGUMENT = (
    "every support atom lies within r bonds of a reacting atom (balls of the predicate radius "
    "around anchors and bounded candidates, which are reacting atoms; embeddings stay in the "
    "ball); reacting atoms of one chain span at most D bonds; closure adds whole pendant "
    "phenyl rings (no backbone extent), radical pi systems (a backbone carbon and its pendant "
    "ring: allylic radicals are outside the domain) and end alkenes (one backbone bond beyond "
    "the support), so reach = r + 1 per side")


def span_bound(constants: FamilyConstants, provider=None) -> SpanBound:
    """The family's span bound.  Without ``provider``: two backbone atoms per
    repeat (the base values).  With one: ``k`` = its
    ``coordinate_backbone_atoms``, ``repeats = ceil((A + k - 1)/k)`` and ``L =
    floor((A + k - 1)/k)``."""
    reach = int(constants.radius) + CLOSURE_EXTENSION
    atoms = int(constants.reacting_distance) + 2 * reach + 1
    if provider is None:
        repeats = math.ceil((atoms + 1) / 2)
        L = (atoms + 1) // 2
        formula = (f"reach = r + {CLOSURE_EXTENSION} = {reach}; backbone atoms A = D + 2*reach "
                   f"+ 1 = {constants.reacting_distance} + {2 * reach} + 1 = {atoms}; repeats = "
                   f"ceil((A + 1)/2) = {repeats}; L = floor((A + 1)/2) = {L} (largest n with "
                   f"2n - 1 <= A: one end on each repeat phase; in-phase chains need 2n <= A)")
    else:
        provider = provider_object(provider)
        k = provider.coordinate_backbone_atoms
        if not k:
            raise WindowError(f"refused:span: provider {provider.name!r} has no backbone")
        repeats = -(-(atoms + k - 1) // k)
        L = (atoms + k - 1) // k
        formula = (f"reach = r + {CLOSURE_EXTENSION} = {reach}; backbone atoms A = D + 2*reach "
                   f"+ 1 = {constants.reacting_distance} + {2 * reach} + 1 = {atoms}; k = "
                   f"{k} backbone atoms per {provider.name} coordinate unit; repeats = "
                   f"ceil((A + k - 1)/k) = {repeats}; L = floor((A + k - 1)/k) = {L}")
    return SpanBound(constants.family, int(constants.radius), int(constants.reacting_distance),
                     reach, atoms, repeats, L, formula, SPAN_ARGUMENT)


def chain_backbone_atoms(spec) -> int:
    """Backbone carbons of a PS chain spec (any phase)."""
    if getattr(spec, "monomer", False):
        return 2
    return 2 * spec.units - int(bool(getattr(spec, "head_shift", False))) + \
        int(bool(getattr(spec, "tail_shift", False)))


def family_constants(family) -> FamilyConstants:
    """``r`` (largest predicate radius over the family's tree nodes) and ``D``
    (largest bond distance between labelled heavy atoms of one connected
    template structure, reactants and products after the recipe), read from
    the loaded RMG family."""
    import networkx as nx
    from rmgpy.data.base import LogicNode
    from rmgpy.kmc import channels as chm
    radius = 0
    for entry in family.groups.entries.values():
        if isinstance(entry.item, LogicNode):
            continue
        for group in chm._group_atoms(family.groups, entry):
            r = chm._group_radius(group)
            if r is not None:
                radius = max(radius, r)
    slots = [chm._group_atoms(family.groups, e) for e in family.forward_template.reactants]
    actions = family.forward_recipe.actions
    distance = 0
    free = False
    for combo in product(*slots):
        graph = nx.Graph()
        labels = defaultdict(list)
        heavy = set()
        for s, group in enumerate(combo):
            index = {a: k for k, a in enumerate(group.atoms)}
            for a in group.atoms:
                node = (s, index[a])
                graph.add_node(node)
                if a.label:
                    labels[a.label].append(node)
                    if not all(t.label == "H" for t in a.atomtype):
                        heavy.add(node)
                for c in a.bonds:
                    graph.add_edge(node, (s, index[c]))
        if any(len(v) > 1 for v in labels.values()):
            free = True

        def ends(l1, l2):
            a, c = labels.get(l1, []), labels.get(l2, [])
            if l1 == l2:
                return [tuple(a[:2])] if len(a) >= 2 else []
            return [(x, y) for x in a for y in c]

        product_graph = graph.copy()
        for action in actions:
            if action[0] == "FORM_BOND":
                for x, y in ends(action[1], action[3]):
                    product_graph.add_edge(x, y)
            elif action[0] == "BREAK_BOND":
                for x, y in ends(action[1], action[3]):
                    if product_graph.has_edge(x, y):
                        product_graph.remove_edge(x, y)
        for g in (graph, product_graph):
            for comp in nx.connected_components(g):
                nodes = [n for n in comp if n in heavy]
                sub = g.subgraph(comp)
                for x, y in combinations(nodes, 2):
                    distance = max(distance, nx.shortest_path_length(sub, x, y))
    directions = ("forward",) + (("reverse",) if family.reverse_template is not None else ())
    molecularity = len(family.forward_template.reactants)
    if family.auto_generated:
        root = slots[0][0] if slots and slots[0] else None
        molecularity = 1 if root is None else nx.number_connected_components(
            _group_graph(root))
    reverse_mol = len(family.forward_template.products)
    bimolecular = molecularity == 2 or (family.reverse_template is not None and reverse_mol == 2)
    dirs = tuple((d, molecularity if d == "forward" else reverse_mol) for d in directions)
    return FamilyConstants(family.label, radius, distance, free, bimolecular, dirs,
                           per_reactant=tuple(d for d, n in dirs
                                              if n == 2 and per_reactant_tree(family, d)),
                           joined=joined_templates(family, dirs))


def _group_graph(group):
    import networkx as nx
    graph = nx.Graph()
    index = {a: k for k, a in enumerate(group.atoms)}
    for a in group.atoms:
        graph.add_node(index[a])
        for c in a.bonds:
            graph.add_edge(index[a], index[c])
    return graph


def window_span(witness, atoms, providers=None) -> dict:
    """Per participant: the window's backbone extent (atoms) and repeat units,
    and the physical ends (port keys) inside the window.  Backbone indices
    come from the participant's provider over its own rows
    (``witness_backbone_indices``)."""
    atoms = set(atoms)
    rows = witness.reactant_rows()
    out = {}
    for participant in witness.participants:
        meta = participant.metadata
        provider = provider_object(meta.provider, providers)
        own = {a: rows[a] for a in participant.atom_ids if a in rows}
        indices = provider.witness_backbone_indices(own, meta)
        idx, units = [], []
        for a, u, _k in meta.coordinates:
            if a in atoms and a in indices:
                idx.append(indices[a])
                units.append(u)
        ends = tuple(sorted(port.port_id for port in meta.ports
                            if not port.artificial and port.atom_id in atoms))
        all_ends = tuple(sorted(port.port_id for port in meta.ports if not port.artificial))
        out[participant.index] = {
            "provider": meta.provider,
            "backbone_atoms": (max(idx) - min(idx) + 1) if idx else 0,
            "units": (max(units) - min(units) + 1) if units else 0,
            "physical_ends_inside": ends, "physical_ends": all_ends,
            "both_ends": len(all_ends) == 2 and len(ends) == 2}
    return out


def repeats_bound(bound: SpanBound, provider) -> int:
    """The bound in the provider's own repeat units (an object or a name):
    ``A`` backbone atoms touch at most ``ceil((A + k - 1) / k)`` repeats of
    ``k = coordinate_backbone_atoms`` backbone atoms (any phase).  A provider
    without a backbone has no bound (refused; ``check_span`` skips it)."""
    provider = provider_object(provider)
    k = provider.coordinate_backbone_atoms
    if not k:
        raise WindowError(f"refused:span: provider {provider.name!r} has no backbone")
    return -(-(bound.backbone_atoms + k - 1) // k)


def check_span(witness, atoms, bound: SpanBound, providers=None) -> dict:
    """``window_span`` with the bound verdict; raises ``SpanBoundError`` when
    a participant's window exceeds the bound: backbone atoms against ``A``,
    repeat units against ``repeats_bound`` in that participant's provider
    units."""
    span = window_span(witness, atoms, providers)
    for index, s in span.items():
        provider = provider_object(s["provider"], providers)
        if not provider.coordinate_backbone_atoms:
            continue    # no backbone (a small species): nothing to bound
        limit = repeats_bound(bound, provider)
        if s["backbone_atoms"] > bound.backbone_atoms or s["units"] > limit:
            raise SpanBoundError(f"{bound.family}: participant {index} ({s['provider']}) window "
                                 f"spans {s['backbone_atoms']} backbone atoms / {s['units']} "
                                 f"units > bound {bound.backbone_atoms} / {limit}")
    return span


@dataclass(frozen=True)
class LadderItem:
    label: str
    kinds: tuple
    dp: int
    specs: tuple
    control: bool = False
    long_index: Optional[int] = None   # the artificial (long) participant of a pair
    dps: tuple = ()                    # the DP of each short participant
    directions: tuple = ()             # pairs: the directions their kinds are run in
    pair_kinds: tuple = ()             # pairs: the participant kinds (sorted)
    classes: tuple = ()                # generated items: the site classes they represent
    leaf_tuples: tuple = ()            # per-reactant pairs: the leaf tuples they must cover
    sites: tuple = ()                  # class pairs: ((spec repr, slot, site), ...) expanded


@dataclass(frozen=True)
class LadderPlan:
    family: str
    L: int
    items: tuple
    excluded: tuple
    pair_kinds: tuple = ()
    single_directions: tuple = ()
    class_cover: tuple = ()   # ((class label, how), ...): classes covered by their closing

    @property
    def excluded_keys(self) -> set:
        return {key for key, _reason in self.excluded}


LONG_UNITS = 4
PARTICIPANT_KINDS = ("alkene", "closed", "radical")


def kind_specs(kind: str, dp: int, provider=None) -> tuple:
    """Representative chains of one participant kind at DP ``dp`` (physical
    ends): ``radical`` a benzylic or a primary end radical, ``closed`` the
    saturated chain, ``alkene`` a vinylidene or a styryl end; at DP 1 also the
    one-carbon species of the kind (benzyl and methyl radicals; toluene and
    methane).  A declared ``provider`` is refused (``_roster_guard``)."""
    _roster_guard(provider)
    S = b.PsChainSpec
    if kind == "radical":
        out = [S(dp, head_end="methyl", tail_end="benzylic_radical"),
               S(dp, head_end="primary_radical", tail_end="benzyl")]
        if dp == 1:
            out += [S(1, head_shift=True, head_end="benzyl", tail_end="benzylic_radical"),
                    S(0, tail_shift=True, head_end="primary_radical", tail_end="methyl")]
    elif kind == "closed":
        out = [S(dp, head_end="methyl", tail_end="benzyl")]
        if dp == 1:
            out += [S(1, head_shift=True, head_end="benzyl", tail_end="benzyl"),
                    S(0, tail_shift=True, head_end="methyl", tail_end="methyl")]
    elif kind == "alkene":
        out = [S(dp, head_end="vinylidene", tail_end="benzyl"),
               S(dp, head_end="methyl", tail_end="styryl")]
    else:
        raise WindowError(f"unknown participant kind {kind!r}")
    return tuple(out)


# the long artificial chains of a kind: every interior reacting-site kind with
# no declared feature within the span bound (both sides continue; contexts
# build the ends at near ports): interior CH2 / CH(Ph) on the closed chain,
# both interior radical kinds (``INTERIOR_RADICALS``: benzylic, secondary),
# an alkene end
LONG_SPECS = {"closed": lambda: (b.PsChainSpec(LONG_UNITS),),
              "radical": lambda: (b.PsChainSpec(LONG_UNITS, radical_units=(LONG_UNITS // 2,)),
                                  b.PsChainSpec(LONG_UNITS,
                                                secondary_radical_units=(LONG_UNITS // 2,))),
              "alkene": lambda: (b.PsChainSpec(LONG_UNITS, head_end="vinylidene"),)}
# each entry refuses a declared provider (``_roster_guard``)
LONG_SPECS = {kind: _roster_entry(make) for kind, make in LONG_SPECS.items()}


def pair_roster(certifier, family, probe_dps=(2, 3), provider=None) -> tuple:
    """The participant-kind combinations a bimolecular direction needs, from
    the family's own reactant templates: RMG applies them (``member_cases``,
    both slot orders) to the representative chains of every kind pair at the
    probe DPs; a pair is in the roster of that direction iff some member has
    every species in the declared domain.  ((direction, (kind, kind)), ...),
    kinds sorted.  A declared ``provider`` is refused."""
    _roster_guard(provider)
    from itertools import combinations_with_replacement
    constants = family_constants(family)
    P = b.PROVIDERS["polystyrene"]
    out = []
    for direction, molecularity in constants.directions:
        if molecularity != 2:
            continue
        for kinds in combinations_with_replacement(PARTICIPANT_KINDS, 2):
            found = False
            for dp in probe_dps:
                for first in kind_specs(kinds[0], dp)[:2]:
                    for second in kind_specs(kinds[1], dp)[:2]:
                        builds = [P.build("k0", first), P.build("k1", second)]
                        if not all(dom.check(x.molecule).accepted for x in builds):
                            continue
                        cases = member_cases(certifier, family, builds, direction)
                        if any(not species_domain_refusal(c.witness) for c in cases):
                            found = True
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                out.append((direction, tuple(kinds)))
    return tuple(out)


def _ends_label(spec):
    if isinstance(spec, tuple):
        return f"{spec[0]}:{spec[1]}"
    head = spec.head_end or "cut"
    tail = spec.tail_end or "cut"
    extra = []
    if spec.head_shift:
        extra.append("hs")
    if spec.tail_shift:
        extra.append("ts")
    for name, flag in (("bz", spec.radical_units), ("sec", spec.secondary_radical_units),
                       ("pl", spec.pendant_loss_units)):
        if flag:
            extra.append(f"{name}{','.join(map(str, flag))}")
    return f"DP{spec.units}[{head}|{tail}{'|' + ','.join(extra) if extra else ''}]"


# ---------------------------------------------------------------------------
# generated roster: reacting-site environment classes (rework 2)
# ---------------------------------------------------------------------------
# Every declared feature of ``domain`` mapped to the site kinds it carries
# (``declared_site_kinds``).  A declared key missing from these maps fails
# the derivation: a new declared feature is never silently left out.
END_SITE_KINDS = {
    "methyl": ("end_methyl",),
    "benzyl": ("end_benzyl",),
    "primary_radical": ("terminal_primary",),
    "benzylic_radical": ("terminal_benzylic",),
    "vinylidene": ("vinylidene_terminal", "vinylidene_partner"),
    "styryl": ("styryl_terminal", "styryl_partner"),
}
RADICAL_SITE_KINDS = {
    "benzylic": ("terminal_benzylic", "interior_benzylic", "benzyl_radical"),
    "primary_end": ("terminal_primary",),
    "secondary_backbone": ("interior_secondary",),
    "methyl": ("methyl_radical",),
    "h_atom": ("h_atom",),
}
# radical kinds that can sit on an interior backbone carbon, and the spec
# field that places them (benzylic: unit position 1; secondary: position 0)
INTERIOR_RADICALS = {"benzylic": "radical_units", "secondary_backbone": "secondary_radical_units"}
SMALL_SITE_KINDS = {"H": ("h_atom",), "H2": ("h2",)}
# the interior closed-shell carbons of the declared chain (``domain``: a
# saturated backbone, phenyl pendants, every other substituent H)
BACKBONE_SITE_KINDS = ("backbone_ch2", "backbone_chph")
# the one-carbon species ``domain._check_form`` accepts (a lone backbone
# carbon: ``methane``, ``methyl_radical``, ``benzyl``, ``benzylic_radical``)
ONE_CARBON = {
    "methane": lambda: b.PsChainSpec(0, tail_shift=True, head_end="methyl", tail_end="methyl"),
    "methyl_radical": lambda: b.PsChainSpec(0, tail_shift=True, head_end="primary_radical",
                                            tail_end="methyl"),
    "benzyl": lambda: b.PsChainSpec(1, head_shift=True, head_end="benzyl", tail_end="benzyl"),
    "benzylic_radical": lambda: b.PsChainSpec(1, head_shift=True, head_end="benzyl",
                                              tail_end="benzylic_radical"),
}
ONE_CARBON_SITE_KINDS = {"methane": ("methane",), "methyl_radical": ("methyl_radical",),
                         "benzyl": ("toluene_methyl",), "benzylic_radical": ("benzyl_radical",)}
# probe partner candidates of a bimolecular direction, cheapest first (no
# resonance before the phenyl species): a template slot holds the same
# reactant whatever fills the other slot, and the partner's own products do
# not depend on the pattern, so the cheapest declared species that fills a
# slot with in-domain products (``probe_partners``) probes every pattern in
# the other slot
PROBE_CANDIDATES = (
    ("H", lambda: ("small", "H")),
    ("H2", lambda: ("small", "H2")),
    ("methyl_radical", ONE_CARBON["methyl_radical"]),
    ("methane", ONE_CARBON["methane"]),
    ("2-phenylethyl", lambda: b.PsChainSpec(1, head_end="primary_radical", tail_end="benzyl")),
    ("ethylbenzene", lambda: b.PsChainSpec(1, head_end="methyl", tail_end="benzyl")),
    ("styrene", lambda: b.PsChainSpec(1, head_end="vinylidene", tail_end="styryl")),
    ("1-phenylethyl", lambda: b.PsChainSpec(1, head_end="methyl", tail_end="benzylic_radical")),
)
PATTERN_DP_LIMIT = 8
# closing partners (``closing_partners``): the probe candidates, then DP2
# chains of every kind (R_Addition_MultipleBond's one-carbon and DP1 radicals
# all end refused against styrene: an element-count guard, a depository match
# or an unknown depository exclusion)
CLOSING_CANDIDATES = PROBE_CANDIDATES + (
    ("DP2 primary radical", lambda: b.PsChainSpec(2, head_end="primary_radical",
                                                  tail_end="benzyl")),
    ("DP2 closed", lambda: b.PsChainSpec(2, head_end="methyl", tail_end="benzyl")),
    ("DP2 vinylidene", lambda: b.PsChainSpec(2, head_end="vinylidene", tail_end="benzyl")),
    ("DP2 benzylic radical", lambda: b.PsChainSpec(2, head_end="methyl",
                                                   tail_end="benzylic_radical")),
)


def declared_site_kinds() -> dict:
    """{site kind: the ``domain`` declaration that carries it}: every end group
    (``END_GROUPS``), radical kind (``RADICAL_KINDS``), small species
    (``SMALL_SPECIES``), the backbone's closed-shell carbons and the
    one-carbon species.  Raises ``WindowError`` for a declared key without a
    mapping (a new declared feature cannot be silently missed)."""
    out = {}
    for decl, table, keys in (("END_GROUPS", END_SITE_KINDS, dom.END_GROUPS),
                              ("RADICAL_KINDS", RADICAL_SITE_KINDS, dom.RADICAL_KINDS),
                              ("SMALL_SPECIES", SMALL_SITE_KINDS, dom.SMALL_SPECIES)):
        for key in keys:
            if key not in table:
                raise WindowError(f"domain.{decl} declares {key!r}, which the generated roster "
                                  "does not map to a reacting-site kind")
            for kind in table[key]:
                out.setdefault(kind, f"{decl}[{key}]")
    for key in INTERIOR_RADICALS:
        if key not in dom.RADICAL_KINDS:
            raise WindowError(f"interior radical {key!r} is not a declared radical kind")
    for kind in BACKBONE_SITE_KINDS:
        out.setdefault(kind, "backbone")
    for name, kinds in ONE_CARBON_SITE_KINDS.items():
        for kind in kinds:
            out.setdefault(kind, f"one_carbon[{name}]")
    return out


RADICAL_SITES = frozenset(k for v in RADICAL_SITE_KINDS.values() for k in v)


@dataclass(frozen=True, order=True)
class Pattern:
    """A declared participant pattern: ``chain`` (a PS chain: its two end
    groups, sorted, and the kinds of its interior radical sites, sorted),
    ``one_carbon`` or ``small`` (``ends`` holds the species name)."""

    kind: str
    ends: tuple = ()
    interior: tuple = ()

    @property
    def label(self) -> str:
        if self.kind != "chain":
            return f"{self.kind}:{self.ends[0]}"
        extra = f"+{'+'.join(self.interior)}" if self.interior else ""
        return f"{self.ends[0]}|{self.ends[1]}{extra}"

    @property
    def radical_sites(self) -> int:
        if self.kind == "small":
            return int(self.ends[0] == "H")
        if self.kind == "one_carbon":
            return int(any(k in RADICAL_SITES for k in ONE_CARBON_SITE_KINDS[self.ends[0]]))
        return sum(_end_is_radical(e) for e in self.ends) + len(self.interior)

    @property
    def features(self) -> tuple:
        """The declared features of the pattern (end groups, interior radical
        kinds, the species name of a one-carbon or small species)."""
        return tuple(self.ends) + tuple(f"interior:{k}" for k in self.interior)


def _end_is_radical(end) -> bool:
    return any(k in RADICAL_SITES for k in END_SITE_KINDS[end])


def participant_patterns(provider=None) -> tuple:
    """Every declared participant pattern: chains with any two declared end
    groups and up to ``MAX_RADICAL_SITES`` radical sites in total (radical ends
    plus interior radicals of each interior kind), the one-carbon species and
    the small species.  A declared ``provider`` is refused."""
    _roster_guard(provider)
    declared_site_kinds()       # every declared key is mapped
    ends = sorted(dom.END_GROUPS)
    out = []
    for pair in combinations_with_replacement(ends, 2):
        room = dom.MAX_RADICAL_SITES - sum(_end_is_radical(e) for e in pair)
        for m in range(max(room, -1) + 1):
            for interior in combinations_with_replacement(sorted(INTERIOR_RADICALS), m):
                out.append(Pattern("chain", tuple(pair), tuple(interior)))
    out += [Pattern("one_carbon", (name,)) for name in ONE_CARBON]
    out += [Pattern("small", (name,)) for name in dom.SMALL_SPECIES]
    return tuple(out)


def small_build(name: str):
    """A small declared species (``H`` atom, ``H2``) as a provider-free build:
    no ports, no repeat coordinates (no heavy atom)."""
    from rmgpy.kmc import rmg_moves as rm
    if name == "H":
        rows = {0: ("H", None, 0, 1, 0, ())}
    elif name == "H2":
        rows = {0: ("H", None, 0, 0, 0, ((1, 1.0),)), 1: ("H", None, 0, 0, 0, ((0, 1.0),))}
    else:
        raise WindowError(f"undeclared small species {name!r}")
    metadata = b.BoundaryMetadata("small", (), ())
    metadata.check(rows, name)
    return b.ProxyBuild(name, rm.molecule_from_rows(rows), metadata, tuple(sorted(rows)))


def _chain_spec(head, tail, n, **interior):
    return b.PsChainSpec(n, head_end=head, tail_end=tail, head_shift=head in b.TAIL_ENDS,
                         tail_shift=tail in b.HEAD_ENDS, **interior)


def _realizes(pattern, spec) -> bool:
    """The spec builds an in-domain chain whose end groups and radical kinds
    are the pattern's (``domain.check``'s own reading).  A declared spec is
    refused."""
    _roster_guard(spec)
    try:
        build = b.PROVIDERS["polystyrene"].build("pattern", spec)
    except b.BoundaryError:
        return False
    if chain_backbone_atoms(spec) < 2:
        return False
    verdict = dom.check(build.molecule)
    if not verdict.accepted or verdict.kind != "ps_chain":
        return False
    if sorted(verdict.details["ends"]) != sorted(pattern.ends):
        return False
    kinds = sorted(k for _a, k in verdict.details["radical_sites"])
    want = sorted([{"primary_radical": "primary_end", "benzylic_radical": "benzylic"}[e]
                   for e in pattern.ends if _end_is_radical(e)] + list(pattern.interior))
    return kinds == want


def _backbone_position(spec, unit, position) -> int:
    return 2 * unit + position - int(bool(spec.head_shift))


_PATTERN_SPECS = {}


def pattern_specs(pattern, dp: int, provider=None) -> tuple:
    """The canonical specs realizing ``pattern`` at DP ``dp`` (constructible and
    accepted by ``domain.check`` with the pattern's ends and radical kinds).
    A chain's first end group sits at the head; interior radicals take the
    in-domain placements nearest the head and nearest the tail (two radicals:
    the first as near the head, the second as near the tail as the domain
    allows, both kind orders), so the free stretch beside each feature is as
    long as the DP allows.  A declared ``provider`` is refused."""
    _roster_guard(provider)
    key = (pattern, dp)
    if key not in _PATTERN_SPECS:
        _PATTERN_SPECS[key] = _pattern_specs(pattern, dp)
    return _PATTERN_SPECS[key]


def _pattern_specs(pattern, dp: int) -> tuple:
    if pattern.kind == "small":
        return (("small", pattern.ends[0]),) if dp == 1 else ()
    if pattern.kind == "one_carbon":
        return (ONE_CARBON[pattern.ends[0]](),) if dp == 1 else ()
    head, tail = pattern.ends
    if not pattern.interior:
        spec = _chain_spec(head, tail, dp)
        return (spec,) if _realizes(pattern, spec) else ()
    sites = []
    for kinds in sorted(set(permutations(pattern.interior))):
        for units in product(range(dp + 1), repeat=len(kinds)):
            fields = defaultdict(tuple)
            for kind, u in zip(kinds, units):
                fields[INTERIOR_RADICALS[kind]] += (u,)
            spec = _chain_spec(head, tail, dp, **{k: tuple(sorted(v)) for k, v in fields.items()})
            pos = [_backbone_position(spec, u, 1 if kind == "benzylic" else 0)
                   for kind, u in zip(kinds, units)]
            if pos != sorted(set(pos)) or any(y - x <= dom.RADICAL_SEPARATION
                                              for x, y in zip(pos, pos[1:])):
                continue    # the domain refuses such a pair (checked again by _realizes)
            if _realizes(pattern, spec):
                sites.append((tuple(pos), kinds, spec))
    if not sites:
        return ()
    out = []
    if len(pattern.interior) == 1:
        sites.sort(key=lambda s: s[0])
        picks = [sites[0], sites[-1]]
    else:
        picks = []
        for kinds in sorted({s[1] for s in sites}):
            mine = [s for s in sites if s[1] == kinds]
            picks.append(min(mine, key=lambda s: (s[0][0], -s[0][-1])))
    for _pos, _kinds, spec in picks:
        if spec not in out:
            out.append(spec)
    return tuple(out)


def pattern_min_dp(pattern, limit: int = PATTERN_DP_LIMIT) -> Optional[int]:
    """The smallest DP at which the pattern is constructible and in domain."""
    for dp in range(1, limit + 1):
        if pattern_specs(pattern, dp):
            return dp
    return None


def build_spec(spec, name="p"):
    """A provider build of a spec (``PsChainSpec``, ("alkane", n[, head end,
    tail end]) or ("small", name))."""
    # every spec is a witness spec first (refused otherwise); a declared spec
    # builds through its registered provider
    spec = witness_spec(spec)
    if isinstance(spec, DeclaredSpec):
        return b.PROVIDERS[spec.provider].build(name, spec.chain)
    if isinstance(spec, tuple) and spec[0] == "alkane":
        ends = tuple(spec[2:4]) + (None,) * (4 - len(spec))
        return b.PROVIDERS["linear_alkane"].build(spec[1], head_end=ends[0], tail_end=ends[1])
    if isinstance(spec, tuple) and spec[0] == "small":
        return small_build(spec[1])
    return b.PROVIDERS["polystyrene"].build(name, spec)


def site_kind(rows, atom, cyclic=None) -> str:
    """The reacting-site kind of a heavy atom (or a hydrogen-only species) of
    one participant's rows: radical sites (``RADICAL_SITE_KINDS``), alkene end
    carbons, closed-shell carbons (ends, backbone, one-carbon species) and
    ``ring``.  A hydrogen bonded to a heavy atom takes that atom's kind."""
    row = rows[atom]
    if row[0] == "H":
        heavy = [n for n, _o in row[5] if rows[n][0] != "H"]
        if heavy:
            return site_kind(rows, heavy[0], cyclic)
        return "h_atom" if row[3] else "h2"
    cyclic = _cyclic_atoms(rows) if cyclic is None else cyclic
    if atom in cyclic:
        return "ring"
    backbone = [n for n, _o in row[5] if rows[n][0] != "H" and n not in cyclic]
    pendant = any(n in cyclic for n, _o in row[5])
    double = any(float(o) == 2.0 for n, o in row[5] if n in backbone)
    if row[3]:
        if not backbone:
            return "benzyl_radical" if pendant else "methyl_radical"
        if len(backbone) == 1:
            return "terminal_benzylic" if pendant else "terminal_primary"
        return "interior_benzylic" if pendant else "interior_secondary"
    if double:
        if len(backbone) == 1:
            return "styryl_terminal" if pendant else "vinylidene_terminal"
        return "vinylidene_partner" if pendant else "styryl_partner"
    if not backbone:
        return "toluene_methyl" if pendant else "methane"
    if len(backbone) == 1:
        return "end_benzyl" if pendant else "end_methyl"
    return "backbone_chph" if pendant else "backbone_ch2"


def pattern_site_kinds(pattern, spec) -> frozenset:
    """Site kinds present on a pattern's spec (ring carbons excluded); a
    declared spec is refused."""
    _roster_guard(spec)
    from rmgpy.kmc import rmg_moves as rm
    rows = rm.molecule_rows(build_spec(spec).molecule)
    cyclic = _cyclic_atoms(rows)
    kinds = {site_kind(rows, a, cyclic) for a in rows
             if rows[a][0] != "H" or pattern.kind == "small"}
    return frozenset(kinds - {"ring"})


@dataclass(frozen=True)
class SiteClass:
    """A reacting-site environment class: the reacting site (``site``: (label,
    site kind) of the slot's labelled atoms on the participant; RMG-free:
    one site kind) on a participant ``pattern`` (its declared features), in
    template ``slot`` of ``direction`` (RMG-free: ``""`` and -1), with the
    smallest DP at which it is constructible and in domain and the spec
    realizing it there."""

    direction: str
    slot: int
    site: tuple
    pattern: Pattern
    dp: int
    spec: Any = field(compare=False)
    partner: Any = field(default=None, compare=False)

    @property
    def key(self) -> tuple:
        return (self.direction, self.slot, self.site, self.pattern)


def site_classes(limit: int = PATTERN_DP_LIMIT, provider=None) -> tuple:
    """The RMG-free class table: every (site kind, pattern) carried by a
    declared pattern at its smallest constructible in-domain DP.  A declared
    ``provider`` is refused."""
    _roster_guard(provider)
    out = []
    for pattern in participant_patterns():
        dp = pattern_min_dp(pattern, limit)
        if dp is None:
            continue
        for spec in pattern_specs(pattern, dp):
            for kind in sorted(pattern_site_kinds(pattern, spec)):
                cls = SiteClass("", -1, (kind,), pattern, dp, spec)
                if cls.key not in {c.key for c in out}:
                    out.append(cls)
    return tuple(out)


def probe_extra(constants: FamilyConstants) -> int:
    """Repeats probed above a pattern's smallest DP: the domain's distance
    constraints are lower bounds of at most ``RADICAL_SEPARATION + 1`` heavy
    bonds between a feature and the event's atoms (which span at most ``D``
    bonds), so ``ceil((S + 1 + D) / 2)`` more repeats free every placement."""
    return -(-(dom.RADICAL_SEPARATION + 1 + int(constants.reacting_distance)) // 2)


def _slot_sites(case, participant: int) -> dict:
    """{slot: ((label, site kind), ...)} of the labelled atoms on
    ``participant`` over every mapping of the member."""
    owned = set(case.witness.participants[participant].atom_ids)
    rows = {a: r for a, r in case.witness.reactant_rows().items() if a in owned}
    cyclic = _cyclic_atoms(rows)
    out = defaultdict(set)
    for move in case.moves:
        per = defaultdict(list)
        for slot, label, atom in move.labels:
            if atom in owned:
                heavy = rows[atom][0] != "H" or not any(rows[n][0] != "H"
                                                        for n, _o in rows[atom][5])
                if heavy:
                    per[slot].append((label, site_kind(rows, atom, cyclic)))
        for slot, sites in per.items():
            out[slot].add(tuple(sorted(sites)))
    return out


def _heavy_key(part) -> tuple:
    """The structure of one species' rows up to its hydrogen ids (heavy atoms
    keep their ids; H-only species by size and radicals)."""
    heavy = {a for a, row in part.items() if row[0] != "H"}
    if not heavy:
        return ("H-only", len(part), sum(row[3] for row in part.values()))
    return tuple(sorted((a, part[a][:5], sum(1 for n, _o in part[a][5] if n not in heavy),
                         tuple(sorted((n, float(o)) for n, o in part[a][5] if n in heavy)))
                        for a in heavy))


_DOMAIN_CACHE = {}


def _part_in_domain(part) -> bool:
    """``domain.check`` of one species' rows, cached by structure; refused
    first, without the resonance search, when a form-independent condition
    already fails: more unpaired electrons than ``MAX_RADICAL_SITES`` (an
    electron count no resonance form changes) or the separated-radical policy
    (``domain.spin_verdict``, which reads pi systems and so every form
    alike)."""
    from rmgpy.kmc import rmg_moves as rm
    key = _heavy_key(part)
    if key not in _DOMAIN_CACHE:
        if len(_DOMAIN_CACHE) > 20000:
            _DOMAIN_CACHE.clear()
        electrons = sum(int(row[3]) for row in part.values())
        if electrons > dom.MAX_RADICAL_SITES:
            _DOMAIN_CACHE[key] = False
        else:
            molecule = rm.molecule_from_rows(part)
            _DOMAIN_CACHE[key] = (electrons < 2 or dom.spin_verdict(molecule).supported) and \
                dom.check(molecule).accepted
    return _DOMAIN_CACHE[key]


def probe_slots(certifier, family, specs, direction) -> dict:
    """{slot: {site, ...}} participant 0 of ``specs`` fills in in-domain moves
    of ``family`` (``direction``): RMG's own moves (``rmg_moves``, the
    certifier's resonance, the ring-edit filter of the witness's policy as in
    ``member_cases``), every reactant and product species accepted by
    ``domain.check``; sites typed by ``site_kind`` on the built reactant.
    Moves, not channels: the probe needs the slot and the domain verdict of
    each mapping only."""
    from rmgpy.kmc import certificates as cm
    from rmgpy.kmc import channels as chm
    from rmgpy.kmc import rmg_moves as rm
    try:
        molecules, participants = b.witness_reactants(builds_of(specs))
    except b.BoundaryError:
        return {}
    if not all(dom.check(m).accepted for m in molecules):
        return {}
    roots = [cm.aromatic_root(m) for m in molecules] if direction == "forward" else \
        list(molecules)
    policy = _builds_policy(None, participants=participants)
    moves = rm.generate_moves(family, roots, direction, resonance=certifier.resonance,
                              select=ring_change_select(family, roots)
                              if policy.ring_edits == "exclude" else None)
    owned = set(participants[0].atom_ids)
    rows = rm.molecule_rows(molecules[0])
    cyclic = _cyclic_atoms(rows)
    kinds = {}
    out = defaultdict(set)
    for move in moves:
        if move.status not in chm.ALLOWED_STATUSES or move.template is None or \
                move.products is None:
            continue
        if not all(_part_in_domain(p) for p in rm.split_rows(dict(move.products))):
            continue
        per = defaultdict(list)
        for slot, label, atom in move.labels:
            if atom not in owned:
                continue
            heavy = rows[atom][0] != "H" or not any(rows[n][0] != "H" for n, _o in rows[atom][5])
            if heavy:
                if atom not in kinds:
                    kinds[atom] = site_kind(rows, atom, cyclic)
                per[slot].append((label, kinds[atom]))
        for slot, sites in per.items():
            out[slot].add(tuple(sorted(sites)))
    return out


def _in_domain_slots(certifier, family, specs, direction) -> set:
    """Slots participant 0 fills in an in-domain move of ``specs``."""
    return set(probe_slots(certifier, family, specs, direction))


def probe_partners(certifier, family, constants: FamilyConstants) -> dict:
    """{bimolecular direction: {slot: (name, spec)}}: per slot the cheapest
    ``PROBE_CANDIDATES`` species that fills it in an in-domain member with
    some other candidate.  A slot no candidate fills is absent (reported by
    ``probe_pattern`` callers)."""
    out = {}
    for direction, molecularity in constants.directions:
        if molecularity != 2:
            continue
        fills = {}
        for name, make in PROBE_CANDIDATES:
            for _other, make2 in PROBE_CANDIDATES:
                for slot in _in_domain_slots(certifier, family, [make(), make2()], direction):
                    fills.setdefault(slot, (name, make()))
        out[direction] = fills
    return out


def probe_pattern(certifier, family, pattern, constants: FamilyConstants,
                  bound: Optional[SpanBound] = None, partners=None) -> tuple:
    """The site classes ``family``'s templates accept for one pattern, per
    direction and slot, with the smallest DP at which RMG generates an
    in-domain member carrying them (every species of the member accepted by
    ``domain.check``): DPs from the pattern's smallest one up to
    ``probe_extra`` more (at most ``L + 1``), its canonical specs, against the
    partner filling each other slot (``probe_partners``) in a bimolecular
    direction."""
    bound = bound or span_bound(constants)
    partners = probe_partners(certifier, family, constants) if partners is None else partners
    start = pattern_min_dp(pattern)
    if start is None:
        return ()
    top = min(bound.L + 1, start + probe_extra(constants)) if pattern.kind == "chain" else 1
    found = {}
    for dp in range(start, max(top, start) + 1):
        for spec in pattern_specs(pattern, dp):
            for direction, molecularity in constants.directions:
                fills = partners.get(direction, {})
                chosen = [None] if molecularity == 1 else \
                    sorted({v for v in fills.values()}, key=repr)
                for partner in chosen:
                    specs = [spec] + ([partner[1]] if partner else [])
                    for slot, sites in probe_slots(certifier, family, specs,
                                                   direction).items():
                        for site in sites:
                            cls = SiteClass(direction, slot, site, pattern, dp, spec,
                                            partner[0] if partner else None)
                            if cls.key not in found:
                                found[cls.key] = cls
    return tuple(found.values())


def class_roster(certifier, family, patterns=None, constants=None, provider=None) -> tuple:
    """Every site class the family's templates accept per direction and slot
    (``probe_pattern`` over ``patterns``, default every declared pattern).  A
    declared ``provider`` is refused."""
    _roster_guard(provider)
    constants = constants or family_constants(family)
    partners = probe_partners(certifier, family, constants)
    out = []
    for pattern in (patterns if patterns is not None else participant_patterns()):
        out.extend(probe_pattern(certifier, family, pattern, constants, partners=partners))
    return tuple(sorted(out, key=lambda c: (c.direction, c.slot, repr(c.site), c.pattern)))


def class_summary(classes, constants: FamilyConstants) -> dict:
    """The class table per direction and slot (site, pattern, smallest DP),
    class counts and the generated item counts (``class_items``)."""
    molecularity = dict(constants.directions)
    table = defaultdict(list)
    for c in classes:
        table[f"{c.direction}|slot{c.slot}"].append(
            [" ".join(f"{lab}:{kind}" for lab, kind in c.site), c.pattern.label, c.dp,
             c.partner])
    items = class_items(classes, molecularity)
    per = defaultdict(int)
    for direction, members in items:
        per[direction] += 1
    return {"classes": len(classes),
            "class_counts": {k: len(v) for k, v in sorted(table.items())},
            "class_table": {k: sorted(v) for k, v in sorted(table.items())},
            "class_items": dict(per),
            "class_specs": {k: len({repr(c.spec) for c in classes
                                    if f"{c.direction}|slot{c.slot}" == k})
                            for k in sorted(table)}}


def class_label(c) -> str:
    return (f"{c.direction}|slot{c.slot}|{' '.join(f'{lab}:{kind}' for lab, kind in c.site)}"
            f"|{c.pattern.label}")


def ladder_class_gaps(plan, classes) -> list:
    """Classes no planned item represents in their direction (the generated
    roster must reach the ladder whole)."""
    have = defaultdict(set)
    for item in plan.items:
        for label in item.classes:
            have[label].update(item.directions or ("*",))
    for label, _how in getattr(plan, "class_cover", ()):   # covered by its closing
        have[label].add("*")
    out = []
    for c in classes:
        label = class_label(c)
        if c.direction not in have.get(label, ()) and "*" not in have.get(label, ()):
            out.append(f"class {label} (DP{c.dp}) has no ladder item")
    return out


def class_items(classes, molecularity: dict) -> tuple:
    """(direction, specs) of the generated ladder: per unimolecular direction
    every class representative alone; per bimolecular direction the product of
    the slot-0 and slot-1 class representatives (distinct specs per slot)."""
    by = defaultdict(lambda: defaultdict(dict))
    for c in classes:
        by[c.direction][c.slot].setdefault(repr(c.spec), (c.spec, c.dp))
    out = []
    for direction in sorted(by):
        slots = by[direction]
        if molecularity.get(direction, 2) == 1:
            for spec, dp in slots.get(0, {}).values():
                out.append((direction, ((spec, dp),)))
            continue
        for x in slots.get(0, {}).values():
            for y in slots.get(1, {}).values():
                out.append((direction, (x, y)))
    return tuple(out)


# ---------------------------------------------------------------------------
# class equivalence, per-reactant trees and the pair pre-filter (rework 3)
# ---------------------------------------------------------------------------
def per_reactant_tree(family, direction: str = "forward") -> bool:
    """True iff RMG rates ``direction`` of ``family`` from one tree per
    reactant, read from the loaded family: the tree is not auto-generated (an
    auto-generated family has one merged ``Root`` whose predicates read every
    reactant together) and every reactant group of the forward template is
    its own top-level node (``family.groups.top``).  RMG then descends each
    reactant's tree on the reactant that carries its labels alone and looks
    the rate up by the tuple of leaves.  Forward only: a reverse member's rate
    is the forward rule's, read on the forward template's reactants (the
    member's products), so a reverse pair is never reduced."""
    if direction != "forward" or getattr(family, "auto_generated", True):
        return False
    reactants = [e.label for e in family.forward_template.reactants]
    tops = {e.label for e in family.groups.top}
    return len(reactants) >= 2 and len(set(reactants)) == len(reactants) and \
        set(reactants) <= tops


def joined_templates(family, directions) -> tuple:
    """((direction, unpaired delta), ...) of the bimolecular directions whose
    template products form ONE species (every reactant atom ends in it: the
    product groups have one connected component in all), with the net
    unpaired-electron change of the direction's recipe (``GAIN_RADICAL`` -
    ``LOSE_RADICAL``, times the template atoms carrying the action's label:
    R_Recombination's one ``LOSE_RADICAL *`` acts on both ``*`` atoms; bond
    actions move no unpaired electron)."""
    import networkx as nx
    from rmgpy.kmc import channels as chm

    def first(entry):
        try:
            groups = chm._group_atoms(family.groups, entry)
        except Exception:    # an unexpandable logic node: not proved joined
            return None
        return groups[0] if groups else None
    out = []
    for direction, molecularity in directions:
        if molecularity != 2:
            continue
        template = family.forward_template if direction == "forward" else \
            family.reverse_template
        recipe = family.forward_recipe if direction == "forward" else \
            getattr(family, "reverse_recipe", None)
        if template is None or recipe is None:
            continue
        products = [first(e) for e in template.products]
        if any(g is None for g in products) or \
                sum(nx.number_connected_components(_group_graph(g)) for g in products) != 1:
            continue
        reactants = [first(e) for e in template.reactants]
        counts = defaultdict(int)
        for g in reactants:
            for a in (g.atoms if g is not None else ()):
                if a.label:
                    counts[a.label] += 1
        delta = 0
        for action in recipe.actions:
            if action[0] in ("GAIN_RADICAL", "LOSE_RADICAL"):
                n = int(action[2]) * max(1, counts.get(action[1], 1))
                delta += n if action[0] == "GAIN_RADICAL" else -n
        out.append((direction, delta))
    return tuple(out)


_UNPAIRED = {}


def unpaired_electrons(spec) -> int:
    """Unpaired electrons of a spec's built species (any resonance form)."""
    key = repr(spec)
    if key not in _UNPAIRED:
        _UNPAIRED[key] = sum(int(a.radical_electrons) for a in build_spec(spec).molecule.atoms)
    return _UNPAIRED[key]


def pair_prefilter(constants: FamilyConstants, direction: str, first, second) -> Optional[str]:
    """Why EVERY member of the pair (``first``, ``second``) is refused, proved
    without expanding it, or None.  The one proof used: a one-product template
    (``joined_templates``) puts all unpaired electrons of both reactants, plus
    the recipe's delta, on one species in every mapping and every resonance
    form, and ``domain`` refuses a species with more than
    ``MAX_RADICAL_SITES`` unpaired electrons (more sites than allowed, or a
    ``u2`` site).  Two-product templates (H abstraction, disproportionation)
    split the electrons per participant, where the class closings already
    showed each class in domain, so no pair of them is dropped."""
    delta = dict(constants.joined).get(direction)
    if delta is None:
        return None
    n = unpaired_electrons(first) + unpaired_electrons(second) + delta
    if n > dom.MAX_RADICAL_SITES:
        return (f"prefilter:domain: every member's single product species carries {n} "
                f"unpaired electrons > MAX_RADICAL_SITES = {dom.MAX_RADICAL_SITES}")
    return None


def closing_partners(certifier, family, constants: FamilyConstants) -> dict:
    """{bimolecular direction: {slot: (name, spec)}}: per template slot the
    cheapest ``CLOSING_CANDIDATES`` species that fills it in an in-domain member
    reaching a fixed point (the closing partner of every class of the OTHER
    slot).  The probe partner H / H2 is the cheapest filler but its members
    end ``refused:guard: element_count`` (a failed carbon-count guard that a
    provider extension could flip), so the closing partner is the first
    candidate whose member closes.  A species of the family's training
    reactions is never a closing partner (``training_keys``): a depository
    match reads both species whole, so a class closed against one can be
    refused whole-species (1-phenylethyl + H2 under H_Abstraction) and the
    closings lose its leaf; this also excludes H and H2 now that their
    members close (C-3)."""
    from rmgpy.kmc import channels as chm
    out = {}
    training, species = training_keys(family), chm.SpeciesKeys()
    cands = [(name, spec) for name, spec in ((name, make()) for name, make in CLOSING_CANDIDATES)
             if species(builds_of([spec])[0].molecule) not in training]
    for direction, molecularity in constants.directions:
        if molecularity != 2:
            continue
        fills = {}
        for name, spec in cands:
            for _other, spec2 in cands:
                for slot, case in _closing_cases(certifier, family, [spec2, spec], direction):
                    if slot in fills:
                        continue
                    observer = MemberObserver(certifier, family, direction, case.wanted)
                    if fixed_point(case.witness, observer).reached:
                        fills[slot] = (name, spec)
            if len(fills) == 2:
                break
        out[direction] = fills
    return out


def _closing_cases(certifier, family, specs, direction):
    """(slot participant 1 fills, member) of the in-domain supported members."""
    try:
        cases = member_cases(certifier, family, builds_of(specs), direction)
    except b.BoundaryError:
        return []
    out = []
    for case in cases:
        if case.support != "supported" or species_domain_refusal(case.witness):
            continue
        for slot in _slot_sites(case, 1):
            out.append((slot, case))
    return out


@dataclass(frozen=True)
class ClassClosure:
    """A class closed once against its closing partner: every member whose
    slot carries the class's site on participant 0 (``results``: fixed points
    or refusals), the pathways (leaf tuples) and own leaves of its
    fixed-point members, the member statuses, participant 0's unpaired
    electrons and guard-capped element counts (whole-species reads), and
    whether participant 0 is a species of a training reaction (RMG's
    depository match reads both species whole)."""

    label: str
    direction: str
    slot: int
    partner: Optional[str]
    results: tuple = field(default=(), compare=False)
    pathways: tuple = ()
    leaves: tuple = ()
    statuses: tuple = ()
    unpaired: int = 0
    elements: tuple = ()
    training: bool = False

    @property
    def fixed(self) -> tuple:
        return tuple(r for r in self.results if r.status == "fixed_point")


def own_leaf(pathway, slot: int) -> str:
    """The leaf of one template slot (per-reactant trees: one leaf per slot;
    a merged root: the single leaf)."""
    pathway = tuple(pathway)
    return pathway[slot] if len(pathway) > 1 and slot < len(pathway) else pathway[0]


_GUARD_MINIMA = {}


def guard_minima(family) -> dict:
    """{element: the largest element-count minimum of any group of the
    family's tree}: counts at or above it pass every element guard."""
    if family.label not in _GUARD_MINIMA:
        from rmgpy.data.base import LogicNode
        from rmgpy.kmc import channels as chm
        out = defaultdict(int)
        for entry in family.groups.entries.values():
            if isinstance(entry.item, LogicNode):
                continue
            for group in chm._group_atoms(family.groups, entry):
                for el, n in (getattr(group, "elementCount", None) or {}).items():
                    out[str(el)] = max(out[str(el)], int(n))
        _GUARD_MINIMA[family.label] = dict(out)
    return _GUARD_MINIMA[family.label]


_TRAINING = {}


def training_keys(family) -> frozenset:
    """Species keys (modulo resonance) of every species of the family's
    training reactions, both sides (RMG's depository match is a reaction
    isomorphism in either direction).  Each molecule is keyed on a copy with
    its atoms numbered: ``SpeciesKeys`` first looks a molecule up by its
    id-keyed rows, and the depository's molecules carry no atom ids, so
    unnumbered every molecule of one size collapsed onto the first one keyed
    (17 keys for 3119 H_Abstraction entries; H2 and 1-phenylethyl missing)."""
    if family.label not in _TRAINING:
        from rmgpy.kmc import channels as chm
        species = chm.SpeciesKeys()
        out = set()
        try:
            depository = family.get_training_depository()
        except Exception:     # no training depository: nothing can match
            depository = None
        for entry in (depository.entries.values() if depository else ()):
            for s in list(entry.item.reactants) + list(entry.item.products):
                mol = (s.molecule[0] if hasattr(s, "molecule") else s).copy(deep=True)
                for k, atom in enumerate(mol.atoms):
                    atom.id = k
                out.add(species(mol))
        _TRAINING[family.label] = frozenset(out)
    return _TRAINING[family.label]


def _snapshot(case, row, keep):
    """A picklable ContextResult of a closed member (the observation's viewer
    a module-level partial, reacting atoms transported)."""
    from dataclasses import replace
    from functools import partial
    if not keep or not keep[0].reached:
        return ContextResult((), row["status"], row.get("reason", ""))
    fp = keep[0]
    obs = fp.observation
    viewer = partial(certificate_view, obs.certificate, frozenset(case.wanted),
                     transport_reacting=True)
    # a fresh witness (its cached read-only row views do not pickle)
    return ContextResult((), fp.status, fp.reason, replace(fp.witness), fp.window,
                         Observation(obs.support, viewer, obs.certificate))


def close_classes(certifier, family, classes, partner=None, *, bound=None,
                  providers=None, provider=None) -> tuple:
    """Close every class of ``classes`` (one direction, slot and spec) once
    against ``partner`` ((name, spec) filling the other slot; None in a
    unimolecular direction): RMG's members of the spec (and partner), each
    member whose slot carries a class's site on participant 0 run once
    (``run_member``, contexts and all).  (closures, rows).  A declared
    ``provider`` (or class or partner spec) is refused."""
    _roster_guard(provider, *(c.spec for c in classes), partner[1] if partner else None)
    from rmgpy.kmc import channels as chm
    first = classes[0]
    if any((c.direction, c.slot, repr(c.spec)) != (first.direction, first.slot, repr(first.spec))
           for c in classes):
        raise WindowError("close_classes: one direction, slot and spec per call")
    specs = [first.spec] + ([partner[1]] if partner else [])
    label = "close " + "+".join(_ends_label(x) if isinstance(x, b.PsChainSpec) else str(x)
                                for x in specs)
    cases = member_cases(certifier, family, builds_of(specs), first.direction, label=label)
    got, rows = defaultdict(list), []
    for case in cases:
        sites = _slot_sites(case, 0).get(first.slot, set())
        mine = [c for c in classes if c.site in sites]
        if not mine:
            continue
        keep = []
        row = run_member(certifier, family, case, bound=bound, providers=providers,
                         contexts=True, specs=specs, keep=keep)
        row["closes"] = sorted(class_label(c) for c in mine)
        rows.append(row)
        snap = _snapshot(case, row, keep)
        for c in mine:
            got[class_label(c)].append((case.pathway, snap))
    molecule = builds_of([first.spec])[0].molecule
    unpaired = sum(int(a.radical_electrons) for a in molecule.atoms)
    minima = guard_minima(family)
    counts = defaultdict(int)
    for a in molecule.atoms:
        counts[a.element.symbol] += 1
    elements = tuple(sorted((el, min(n, minima.get(el, 0))) for el, n in counts.items()))
    training = chm.SpeciesKeys()(molecule) in training_keys(family)
    out = []
    for c in classes:
        mine = got.get(class_label(c), [])
        fixed = [p for p, s in mine if s.status == "fixed_point"]
        out.append(ClassClosure(class_label(c), c.direction, c.slot,
                                partner[0] if partner else None,
                                tuple(s for _p, s in mine), tuple(sorted(set(fixed))),
                                tuple(sorted({own_leaf(p, c.slot) for p in fixed})),
                                tuple(sorted({s.status for _p, s in mine})), unpaired, elements,
                                training))
    return tuple(out), rows


def _class_graph(witness, atoms):
    """The window's marked graph for class equivalence across witnesses
    (``_collapsed_graph``'s, made witness-free): a reacting atom is marked
    reacting (its labels are compared by the certificate view under the
    isomorphism, ``transport_reacting``), and a chain's orientation is not a
    mark: a physical end is marked physical and an exterior bond names how
    many artificial / physical ports lie beyond it, not which (the same end
    chemistry at a spec's head or tail is the same environment)."""
    import networkx as nx
    nodes, edges = window_labels(witness, atoms)
    rows = witness.reactant_rows()
    support = b.edit_support(witness)
    kind = {(p.index, port.port_id): (p.index, "artificial" if port.artificial else "physical")
            for p in witness.participants for port in p.metadata.ports}

    def norm(label):
        before, after, role, reacting, physical, exterior = label
        physical = None if physical is None else (physical[0], "physical")
        exterior = tuple(sorted(((o1, o2, tuple(sorted(kind.get(k, k) for k in reach)))
                                 for o1, o2, reach in exterior), key=repr))
        return (before, after, role, reacting, physical, exterior)
    kept = {a for a in nodes if rows[a][0] != "H" or a in support}
    folded = defaultdict(list)
    for a in nodes:
        if a in kept:
            continue
        heavy = [n for n, _o in rows[a][5] if n in kept]
        folded[heavy[0] if heavy else None].append(repr(norm(nodes[a])))
    graph = nx.Graph()
    for a in kept:
        graph.add_node(a, label=repr((norm(nodes[a]), "r" if a in support else "",
                                      tuple(sorted(folded.get(a, ()))))))
    graph.graph["orphans"] = tuple(sorted(folded.get(None, ())))
    for (x, y), label in edges.items():
        if x in kept and y in kept:
            graph.add_edge(x, y, label=repr(label))
    return graph


def class_window_isomorphisms(first, atoms_first, second, atoms_second,
                              limit=ISOMORPHISM_LIMIT):
    """Mark-preserving isomorphisms (second -> first) of two windows on
    different witnesses (``_class_graph``)."""
    from networkx.algorithms import isomorphism
    g1, g2 = _class_graph(first, atoms_first), _class_graph(second, atoms_second)
    if g1.graph["orphans"] != g2.graph["orphans"] or len(g1) != len(g2) or \
            sorted(d["label"] for _a, d in g1.nodes(data=True)) != \
            sorted(d["label"] for _a, d in g2.nodes(data=True)):
        return []
    matcher = isomorphism.GraphMatcher(
        g2, g1, node_match=lambda x, y: x["label"] == y["label"],
        edge_match=lambda x, y: x["label"] == y["label"])
    out = []
    for mapping in matcher.isomorphisms_iter():
        out.append(dict(mapping))
        if len(out) >= limit:
            break
    return out


def ball_window(result, participant: int, reach: int) -> frozenset:
    """The closed window of ``result`` widened by every atom of
    ``participant`` within ``reach`` heavy bonds of its reacting atoms (the
    span bound when ``reach`` is the tree's largest predicate radius ``r``:
    no predicate of any descent, whatever the partner, reads beyond it),
    closed again (rings, pi systems, hydrogens)."""
    witness = result.witness
    rows = witness.reactant_rows()
    owned = set(witness.participants[participant].atom_ids)
    reacting = set(result.window.support.reacting) & owned
    seen, frontier, depth = set(reacting), set(reacting), 0
    while frontier and depth < reach:
        nxt = set()
        for a in frontier:
            for n, _o in rows[a][5]:
                if n in owned and n not in seen and rows[n][0] != "H":
                    nxt.add(n)
        seen |= nxt
        frontier = nxt
        depth += 1
    window = close(witness, Support(frozenset(set(result.window.atoms) | seen),
                                    frozenset(result.window.support.reacting)))
    if window.refusal:
        raise WindowError(f"ball window: {window.refusal}")
    return window.atoms


def class_results_merge(first: ContextResult, second: ContextResult, *, reach=None,
                        participant: int = 0) -> bool:
    """Two class members' closed windows are equal marked graphs (across
    witnesses) and their certificate views are equal under that
    correspondence; with ``reach`` (a merged-root tree) the windows widened to
    the span bound around participant 0's reacting atoms must be equal too
    (``ball_window``)."""
    if first.status != "fixed_point" or second.status != "fixed_point":
        return False
    windows = [(first.window.atoms, second.window.atoms)]
    if reach is not None:
        windows.append((ball_window(first, participant, reach),
                        ball_window(second, participant, reach)))
    for w1, w2 in windows:
        ok = False
        for mapping in class_window_isomorphisms(first.witness, w1, second.witness, w2):
            if first.observation.view(w1, {a: a for a in w1}) == \
                    second.observation.view(w2, mapping):
                ok = True
                break
        if not ok:
            return False
    return True


def classes_merge(first: ClassClosure, second: ClassClosure, *, reach=None) -> bool:
    """Two classes of one direction and slot are the same rate-selection
    environment: closed against the same partner, the same leaves (pathways
    of their fixed-point members, checked explicitly), the same member
    statuses, and every fixed-point member of each merges with one of the
    other (``class_results_merge``).  A class without a fixed-point member
    never merges.  With ``reach`` (merged-root tree): also equal unpaired
    electrons and guard-capped element counts (whole-species reads)."""
    if (first.direction, first.slot, first.partner) != \
            (second.direction, second.slot, second.partner):
        return False
    if first.pathways != second.pathways:
        return False
    if first.statuses != second.statuses:
        return False
    if reach is not None and (first.unpaired, first.elements) != \
            (second.unpaired, second.elements):
        return False
    fixed1, fixed2 = first.fixed, second.fixed
    if not fixed1 or not fixed2:
        return False
    for xs, ys, flip in ((fixed1, fixed2, False), (fixed2, fixed1, True)):
        for x in xs:
            if not any(class_results_merge(y, x, reach=reach) if flip else
                       class_results_merge(x, y, reach=reach) for y in ys):
                return False
    return True


@dataclass(frozen=True)
class ClassMerge:
    """The class -> representative map of a family's roster (``merge_classes``)
    and what the plan needs of each class: its own leaves, whether its
    closing reached a fixed point, whether its participant is a training
    species."""

    representative: tuple = ()      # ((class label, representative label), ...)
    leaves: tuple = ()              # ((class label, (leaf, ...)), ...)
    fixed: frozenset = frozenset()  # class labels whose closing reached a fixed point
    training: frozenset = frozenset()
    closed: frozenset = frozenset()  # class labels with a closing
    statuses: tuple = ()             # ((class label, (member status, ...)), ...)

    def rep(self, label: str) -> str:
        return dict(self.representative).get(label, label)

    def leaves_of(self, label: str) -> tuple:
        return dict(self.leaves).get(label, ())


def _class_cost(c) -> tuple:
    try:
        heavy = sum(1 for a in build_spec(c.spec).molecule.atoms if a.element.symbol != "H")
    except b.BoundaryError:
        heavy = 10 ** 6
    return (c.dp, heavy, repr(c.spec), class_label(c))


def merge_classes(classes, closures, *, merged_root=(), reach=None) -> ClassMerge:
    """Merge the classes of each direction and slot whose closings are the
    same rate-selection environment (``classes_merge``): a class merges into
    the first representative (cheapest first: smallest DP, fewest heavy
    atoms) it equals, else it is a representative.  ``merged_root``: the
    directions whose tree reads every reactant together (compared with the
    windows widened to ``reach``).  A class without a closing is its own
    representative."""
    by = {c.label: c for c in closures}
    reps = defaultdict(list)
    pairs = []
    for c in sorted(classes, key=_class_cost):
        label = class_label(c)
        mine = by.get(label)
        rep = label
        if mine is not None:
            r = reach if c.direction in merged_root else None
            for other in reps[(c.direction, c.slot)]:
                if classes_merge(by[other], mine, reach=r):
                    rep = other
                    break
            else:
                reps[(c.direction, c.slot)].append(label)
        pairs.append((label, rep))
    return ClassMerge(tuple(sorted(pairs)),
                      tuple(sorted((c.label, c.leaves) for c in closures)),
                      frozenset(c.label for c in closures if c.fixed),
                      frozenset(c.label for c in closures if c.training),
                      frozenset(c.label for c in closures),
                      tuple(sorted((c.label, c.statuses) for c in closures)))


def merge_summary(classes, merge) -> dict:
    """Per direction and slot: class and representative counts, the
    representatives' distinct specs, and the map representative -> classes."""
    by_label = {class_label(c): c for c in classes}
    table = defaultdict(lambda: defaultdict(list))
    for c in classes:
        label = class_label(c)
        table[f"{c.direction}|slot{c.slot}"][merge.rep(label) if merge else label].append(label)
    out = {}
    for key, reps in sorted(table.items()):
        out[key] = {"classes": sum(len(v) for v in reps.values()),
                    "representatives": len(reps),
                    "representative_specs": len({repr(by_label[r].spec) for r in reps
                                                 if r in by_label}),
                    "fixed": sum(1 for r in reps if merge and r in merge.fixed),
                    "map": {r: sorted(v) for r, v in sorted(reps.items())}}
    return out


def leaf_cover(slot_specs, prefilter=None) -> tuple:
    """Spec pairs covering every leaf tuple that occurs: ``slot_specs`` is
    ({spec key: leaves} of slot 0, {spec key: leaves} of slot 1); the tuples
    are slot-0 leaves x slot-1 leaves; pairs are picked greedily (most
    uncovered tuples first, then the cheapest key order), skipping pairs
    ``prefilter(k0, k1)`` drops.  (((k0, k1, tuples), ...), uncovered
    tuples)."""
    first, second = slot_specs
    want = {(x, y) for v in first.values() for x in v for w in second.values() for y in w}
    keys0, keys1 = list(first), list(second)
    chosen, left = [], set(want)
    while left:
        best, gain = None, set()
        for k0 in keys0:
            for k1 in keys1:
                got = {(x, y) for x in first[k0] for y in second[k1]} & left
                if len(got) > len(gain) and not (prefilter and prefilter(k0, k1)):
                    best, gain = (k0, k1), got
        if best is None:
            break
        chosen.append((best[0], best[1], tuple(sorted(gain))))
        left -= gain
    return tuple(chosen), tuple(sorted(left))


def site_select(case, specs, sites) -> bool:
    """A class-pair member is expanded iff, in every slot the pair constrains,
    some participant carries one of the planned representative sites
    (``sites``: (spec repr, slot, site), ``LadderItem.sites``); the other
    members of the two species belong to other classes, each expanded with
    its own representative."""
    allowed = defaultdict(set)
    for key, slot, site in sites:
        allowed[int(slot)].add((key, tuple(tuple(x) for x in site)))
    keys = [repr(s) for s in specs]
    for slot, want in allowed.items():
        if not any((key, site) in want for p, key in enumerate(keys)
                   for site in _slot_sites(case, p).get(slot, ())):
            return False
    return True


def leaf_tuple_gaps(plan, rows) -> list:
    """Every leaf tuple a per-reactant pair item was planned for (``leaf_tuples``)
    needs a fixed-point member with that pathway in its direction."""
    have = defaultdict(set)
    for r in rows:
        if r.get("status") == "fixed_point" and r.get("pathway"):
            have[r.get("direction")].add(tuple(r["pathway"]))
    out = []
    for item in plan.items:
        for direction in item.directions:
            for leaves in item.leaf_tuples:
                if tuple(leaves) not in have[direction]:
                    out.append(f"{direction} leaf tuple {';'.join(leaves)} ({item.label}): no "
                               "fixed-point member")
    return out


def closure_gaps(classes, merge, directions=None) -> list:
    """Every class must have a closing with at least one member carrying its
    site (the probe found it there; ``directions``: the closed directions,
    default all).  A class whose closing reached no fixed
    point (every member refused by status) is accounted, not a gap: it is
    its own representative and is expanded unreduced (``class_plan``)."""
    out = []
    if merge is None:
        return out
    members = dict(merge.statuses)
    for c in classes:
        if directions is not None and c.direction not in directions:
            continue
        label = class_label(c)
        if label not in merge.closed:
            out.append(f"class {label}: no closing")
        elif not members.get(label):
            out.append(f"class {label}: its closing has no member carrying the site")
    return out


def ladder_plan(constants: FamilyConstants, bound: Optional[SpanBound] = None,
                provider=None) -> LadderPlan:
    """The short-chain ladder of a family up to L, plus the DP L + 1 controls:
    intact in-phase chains in every declared end chemistry and the one-carbon
    species (DP 1); scission fragments of both backbone cut phases (phase A:
    in-phase fragments with the cut end radical or alkene; phase B: shifted
    fragments); pendant loss and C-H cleavage products (one per interior unit
    position, saturated ends).  For every bimolecular direction, every kind
    combination of its roster (``constants.pair_kinds``, ``pair_roster``):
    short/short pairs at EVERY DP pair (a, b), 1 <= a, b <= L + 1 (any DP
    above L is a control), each kind's representative chains
    (``kind_specs``), and short/long pairs (each short kind against the
    other kind's long artificial chain, ``LONG_SPECS``).  Species outside the
    declared domain are excluded with the domain's reason; species built
    twice, and pairs of the same two species, are merged (kinds and
    directions joined).  A declared ``provider`` is refused."""
    _roster_guard(provider)
    bound = bound or span_bound(constants)
    L = bound.L
    if constants.classes is None:
        raise WindowError(f"{constants.family}: the generated class roster was not computed "
                          "(FamilyConstants.classes: class_roster)")
    if any(n == 2 for _d, n in constants.directions) and not constants.pair_kinds:
        raise WindowError(f"{constants.family}: a bimolecular family without a pair roster "
                          "(FamilyConstants.pair_kinds: pair_roster)")
    items, excluded = {}, []
    by_species = {}
    keys = {}

    def species_key(spec):
        from rmgpy.kmc import channels as chm
        r = repr(spec)
        if r not in keys:
            try:
                build = build_spec(spec, "ladder")
            except b.BoundaryError as error:
                keys[r] = (None, f"not constructible: {error}")
            else:
                verdict = dom.check(build.molecule)
                keys[r] = (chm.SpeciesKeys()(build.molecule), None) if verdict.accepted else \
                    (None, f"outside the declared domain: {verdict.reason}")
        return keys[r]

    def add(key, kind, spec, dp, classes=(), sites=()):
        skey, why = species_key(spec)
        if why:
            excluded.append((key, why))
            return
        if skey in by_species:
            label = by_species[skey]
            item = items[label]
            if kind not in item.kinds or classes:
                # site selection only between class entries; a grid item of
                # the same species runs every member
                both = set(item.kinds) | {kind} == {"class"}
                items[label] = LadderItem(item.label, item.kinds + ((kind,) if kind not in
                                                                    item.kinds else ()),
                                          item.dp, item.specs, item.control, dps=item.dps,
                                          classes=item.classes + tuple(classes),
                                          sites=tuple(sorted(set(item.sites) | set(sites)))
                                          if both else ())
            if key[0] != item.label and not classes:
                excluded.append((key, f"duplicate species of {label}"))
            return
        label = _ends_label(spec)
        by_species[skey] = label
        items[label] = LadderItem(label, (kind,), dp, (spec,), dp > L, dps=(dp,),
                                  classes=tuple(classes), sites=tuple(sorted(sites)))

    S = b.PsChainSpec
    for dp in range(1, L + 2):
        for head in b.HEAD_ENDS:
            for tail in b.TAIL_ENDS:
                spec = S(dp, head_end=head, tail_end=tail)
                add((dp, head, tail), "intact", spec, dp)
        if dp == 1:      # one-carbon chains: a lone position-1 or position-0 carbon
            for head in b.TAIL_ENDS:
                for tail in b.TAIL_ENDS:
                    add(("C1", head, tail), "intact",
                        S(1, head_shift=True, head_end=head, tail_end=tail), 1)
            for head in b.HEAD_ENDS:
                for tail in b.HEAD_ENDS:
                    add(("C0", head, tail), "intact",
                        S(0, tail_shift=True, head_end=head, tail_end=tail), 1)
        for head in b.HEAD_ENDS:     # phase A: a cut between repeats
            for tail in ("benzylic_radical", "styryl"):
                add(("A", dp, head, tail), "scission:phase_A",
                    S(dp, head_end=head, tail_end=tail), dp)
        for head in ("primary_radical", "vinylidene"):
            for tail in b.TAIL_ENDS:
                add(("A", dp, head, tail), "scission:phase_A",
                    S(dp, head_end=head, tail_end=tail), dp)
        for head in b.HEAD_ENDS:     # phase B: a cut inside a repeat
            for tail in b.HEAD_ENDS:
                add(("B", dp, head, tail), "scission:phase_B",
                    S(dp, head_end=head, tail_end=tail, tail_shift=True), dp)
        for head in b.TAIL_ENDS:
            for tail in b.TAIL_ENDS:
                add(("B", dp, head, tail), "scission:phase_B",
                    S(dp, head_end=head, tail_end=tail, head_shift=True), dp)
        for u in range(dp):
            if 2 * u + 1 < 2 * dp - 1:
                add(("pl", dp, u), "pendant_loss",
                    S(dp, head_end="methyl", tail_end="benzyl", pendant_loss_units=(u,)), dp)
                add(("bz", dp, u), "ch_cleavage:benzylic",
                    S(dp, head_end="methyl", tail_end="benzyl", radical_units=(u,)), dp)
            if u >= 1:
                add(("sec", dp, u), "ch_cleavage:secondary",
                    S(dp, head_end="methyl", tail_end="benzyl", secondary_radical_units=(u,)),
                    dp)
    out = list(items.values())
    molecularity = dict(constants.directions)
    roster = tuple((d, tuple(sorted(k))) for d, k in constants.pair_kinds
                   if molecularity.get(d, 2) == 2)
    pairs, long_pairs = {}, {}

    def label_of(kind, spec):
        return f"{kind}:{_ends_label(spec)}"
    for direction, kinds in roster:
        for a in range(1, L + 2):
            for c in range(1, L + 2):
                for first in kind_specs(kinds[0], a):
                    for second in kind_specs(kinds[1], c):
                        k1, why1 = species_key(first)
                        k2, why2 = species_key(second)
                        if why1 or why2:
                            excluded.append((("pair", direction, kinds, a, c, repr(first),
                                              repr(second)), why1 or why2))
                            continue
                        key = tuple(sorted((k1, k2)))
                        item = pairs.get(key)
                        if item is None:
                            pairs[key] = LadderItem(
                                f"pair DP{a}+DP{c} {label_of(kinds[0], first)}+"
                                f"{label_of(kinds[1], second)}", ("pair:short_short",),
                                max(a, c), (first, second), max(a, c) > L, None, (a, c),
                                (direction,), kinds)
                        elif direction not in item.directions:
                            pairs[key] = LadderItem(item.label, item.kinds, item.dp, item.specs,
                                                    item.control, None, item.dps,
                                                    item.directions + (direction,),
                                                    item.pair_kinds)
        for short, long_ in sorted({(kinds[0], kinds[1]), (kinds[1], kinds[0])}):
            for partner in LONG_SPECS[long_]():
                for a in range(1, L + 2):
                    for first in kind_specs(short, a):
                        k1, why = species_key(first)
                        if why or not species_key(partner)[0]:
                            excluded.append((("long", direction, short, long_, a, repr(first)),
                                             why or species_key(partner)[1]))
                            continue
                        key = (k1, long_, species_key(partner)[0])
                        item = long_pairs.get(key)
                        if item is None:
                            long_pairs[key] = LadderItem(
                                f"pair DP{a} {label_of(short, first)}+long "
                                f"{label_of(long_, partner)}",
                                ("pair:short_long",), a, (first, partner), a > L, 1, (a,),
                                (direction,), tuple(sorted((short, long_))))
                        elif direction not in item.directions:
                            long_pairs[key] = LadderItem(item.label, item.kinds, item.dp,
                                                         item.specs, item.control, 1, item.dps,
                                                         item.directions + (direction,),
                                                         item.pair_kinds)
    # the generated roster (``class_plan``): per unimolecular direction every
    # representative class alone (none when the closings ran: they are those
    # runs); per bimolecular direction the product of the slot
    # representatives, or (per-reactant trees) the leaf-tuple cover and the
    # training-species product; merged with the grid's items
    cover = []
    entries, class_excluded = _class_plan(constants)
    for entry in entries:
        direction, members, covered, kind = entry[:4]
        leaf_tuples = entry[4] if len(entry) > 4 else ()
        sites = tuple(sorted(entry[5])) if len(entry) > 5 else ()
        if kind == "cover":
            cover.extend((label, members) for label in covered)
            continue
        if len(members) == 1:
            (spec, dp), = members
            add(("class", direction, repr(spec)), "class", spec, dp, tuple(sorted(covered)),
                sites)
            continue
        (first, a), (second, c) = members
        k1, why1 = species_key(first)
        k2, why2 = species_key(second)
        if why1 or why2:
            excluded.append((("class pair", direction, repr(first), repr(second)), why1 or why2))
            continue
        key = tuple(sorted((k1, k2)))
        item = pairs.get(key)
        if item is None:
            pairs[key] = LadderItem(f"class pair DP{a}+DP{c} {_ends_label(first)}+"
                                    f"{_ends_label(second)}", (kind,), max(a, c),
                                    (first, second), max(a, c) > L, None, (a, c),
                                    (direction,), ("class",), tuple(sorted(covered)),
                                    tuple(leaf_tuples), sites)
        else:
            # site selection only between class-product entries; anything else
            # on the same species pair runs every member
            both = set(item.kinds) | {kind} == {"pair:class"}
            pairs[key] = LadderItem(item.label, item.kinds + ((kind,) if kind not in item.kinds
                                                              else ()),
                                    item.dp, item.specs, item.control, None, item.dps,
                                    item.directions + ((direction,) if direction not in
                                                       item.directions else ()),
                                    item.pair_kinds, item.classes + tuple(sorted(covered)),
                                    item.leaf_tuples + tuple(leaf_tuples),
                                    tuple(sorted(set(item.sites) | set(sites))) if both else ())
    excluded.extend(class_excluded)
    out = list(items.values())
    out += list(pairs.values()) + list(long_pairs.values())
    singles = tuple(d for d, n in constants.directions if n == 1)
    return LadderPlan(constants.family, L, tuple(out), tuple(excluded), roster, singles,
                      tuple(sorted(cover)))


def _spec_name(spec) -> str:
    return _ends_label(spec) if isinstance(spec, b.PsChainSpec) else repr(spec)


def class_plan(constants: FamilyConstants) -> list:
    """The generated ladder entries: (direction, ((spec, dp), ...), covered
    class labels, kind[, leaf tuples]).  ``constants.merge`` (``merge_classes``)
    names each class's representative; without it every class is its own.

    * unimolecular direction: every representative's spec alone (``class``);
      with closings these ARE the runs, so the classes are covered by them
      (kind ``cover``) and nothing is planned twice;
    * bimolecular direction whose tree reads every reactant together (not in
      ``constants.per_reactant``): the product of the slot-0 and slot-1
      representative specs (``pair:class``), less the pairs
      ``pair_prefilter`` proves fully refused (``class_plan_excluded``);
    * bimolecular per-reactant direction with closings: the classes are
      covered by their representatives' closings (``cover``); one pair per
      leaf tuple that occurs (``leaf_cover`` over the representatives'
      leaves: ``pair:leaf``, each item naming the tuples it must reach); and
      the product of the classes whose species is a training species
      (``pair:training``: RMG's depository match reads both species whole)."""
    return _class_plan(constants)[0]


def _class_plan(constants):
    classes = constants.classes or ()
    merge = constants.merge
    molecularity = dict(constants.directions)
    labels = {class_label(c): c for c in classes}

    def rep(c):
        return labels.get(merge.rep(class_label(c)), c) if merge is not None else c
    slots = defaultdict(lambda: defaultdict(dict))     # direction -> slot -> spec key -> (spec, dp)
    covered = defaultdict(set)                         # (direction, slot, spec key) -> labels
    for c in classes:
        r = rep(c)
        slots[c.direction][r.slot].setdefault(repr(r.spec), (r.spec, r.dp))
        covered[(c.direction, r.slot, repr(r.spec))].add(class_label(c))
    entries, excluded = [], []
    for direction in sorted(slots):
        mine = [c for c in classes if c.direction == direction]
        closed = merge is not None and all(class_label(c) in merge.closed for c in mine)
        by = slots[direction]
        if molecularity.get(direction, 2) == 1:
            if closed:
                entries.append((direction, "closing", tuple(class_label(c) for c in mine),
                                "cover"))
                continue
            # a spec runs the members carrying one of its own classes' sites
            # (``site_select``): its other sites are classes of a smaller DP,
            # each run on its own spec
            own = defaultdict(set)
            for c in mine:
                own[repr(c.spec)].add(c.site)
            for key, (spec, dp) in sorted(by.get(0, {}).items()):
                entries.append((direction, ((spec, dp),), covered[(direction, 0, key)], "class",
                                (), {(key, 0, site) for site in own[key]}))
            continue
        if direction in constants.per_reactant and closed:
            entries.append((direction, "closing", tuple(class_label(c) for c in mine), "cover"))
            leaves = ({}, {})
            for slot in (0, 1):
                for c in mine:
                    if c.slot != slot or merge.rep(class_label(c)) != class_label(c) or \
                            class_label(c) not in merge.fixed:
                        continue
                    got = leaves[slot].setdefault(repr(c.spec), set())
                    got.update(merge.leaves_of(class_label(c)))
            specs = {(s, k): v for s in (0, 1) for k, v in by.get(s, {}).items()}

            def drop(k0, k1):
                return pair_prefilter(constants, direction, specs[(0, k0)][0],
                                      specs[(1, k1)][0])
            chosen, left = leaf_cover(leaves, drop)
            for k0, k1, tuples in chosen:
                entries.append((direction, (specs[(0, k0)], specs[(1, k1)]),
                                covered[(direction, 0, k0)] | covered[(direction, 1, k1)],
                                "pair:leaf", tuples))
            for t in left:
                excluded.append((("leaf tuple", direction, t), "no representative pair "
                                 "survives the pre-filter (every member refused)"))
            # a class whose closing reached no fixed point (every member refused
            # by status) has no leaf to cover: it stays unreduced, paired with
            # every representative spec of the other slot
            seen = set()
            for c in mine:
                if merge.rep(class_label(c)) != class_label(c) or class_label(c) in merge.fixed:
                    continue
                other = 1 - c.slot
                for k, y in sorted(by.get(other, {}).items()):
                    pair = ((c.spec, c.dp), y) if c.slot == 0 else (y, (c.spec, c.dp))
                    key = (repr(pair[0][0]), repr(pair[1][0]))
                    if key in seen:
                        continue
                    seen.add(key)
                    why = pair_prefilter(constants, direction, pair[0][0], pair[1][0])
                    if why:
                        excluded.append((("unclosed pair", direction) + key, why))
                        continue
                    entries.append((direction, pair, {class_label(c)} |
                                    covered[(direction, other, k)], "pair:unclosed"))
            train = [[c for c in mine if c.slot == s and class_label(c) in merge.training]
                     for s in (0, 1)]
            seen = set()
            for x in train[0]:
                for y in train[1]:
                    key = (repr(x.spec), repr(y.spec))
                    if key in seen:
                        continue
                    seen.add(key)
                    why = pair_prefilter(constants, direction, x.spec, y.spec)
                    if why:
                        excluded.append((("training pair", direction) + key, why))
                        continue
                    entries.append((direction, ((x.spec, x.dp), (y.spec, y.dp)),
                                    {class_label(x), class_label(y)}, "pair:training"))
            continue
        # the representatives' own sites per slot and spec: a pair expands the
        # members carrying a representative site in each slot (``site_select``)
        rep_sites = defaultdict(set)
        for c in mine:
            if rep(c) is c or class_label(rep(c)) == class_label(c):
                rep_sites[(c.slot, repr(c.spec))].add(c.site)
        for k0, x in sorted(by.get(0, {}).items()):
            for k1, y in sorted(by.get(1, {}).items()):
                why = pair_prefilter(constants, direction, x[0], y[0])
                if why:
                    # accounted: every member refused, proved (``pair_prefilter``)
                    excluded.append((("class pair", direction, k0, k1), why))
                    entries.append((direction, "prefilter", tuple(
                        covered[(direction, 0, k0)] | covered[(direction, 1, k1)]), "cover"))
                    continue
                sites = {(k0, 0, site) for site in rep_sites[(0, k0)]} | \
                    {(k1, 1, site) for site in rep_sites[(1, k1)]}
                entries.append((direction, (x, y),
                                covered[(direction, 0, k0)] | covered[(direction, 1, k1)],
                                "pair:class", (), sites))
    return entries, excluded


def ladder_direction_gaps(plan: LadderPlan, rows) -> list:
    """Coverage per direction: for every bimolecular direction and kind
    combination of the roster, and every unimolecular direction (single
    chains), some member at DP L and at DP L + 1 (the largest short DP of its
    pair) reached a fixed point.  ``rows``: run_member rows with
    ``direction``, ``pair_kinds`` (pairs) and ``dps``."""
    have = set()
    for r in rows:
        if r.get("status") != "fixed_point" or not r.get("dps"):
            continue
        kinds = tuple(r["pair_kinds"]) if r.get("pair_kinds") else None
        have.add((r.get("direction"), kinds, max(d for d in r["dps"] if d is not None)))
    need = [(d, tuple(k)) for d, k in plan.pair_kinds] + [(d, None) for d in
                                                          plan.single_directions]
    gaps = []
    for direction, kinds in need:
        for dp in (plan.L, plan.L + 1):
            if (direction, kinds, dp) not in have:
                gaps.append(f"{direction} {'+'.join(kinds) if kinds else 'single chain'} "
                            f"DP{dp} (L={plan.L}): no member reached a fixed point")
    return gaps


# ---------------------------------------------------------------------------
# spin and resonance diagnostics (item 7)
# ---------------------------------------------------------------------------
def diagnostics(witness, atoms, policy=None) -> dict:
    """Radical sites in the window (both sides), the separated-radical verdict
    (``policy``'s domain, default ``LEGACY_V1``) of every species touched and
    the matcher's resonance form count per whole component used (or the
    overflow cutoff, never reported as a count)."""
    from rmgpy.kmc import rmg_moves as rm
    from rmgpy.kmc import window_matcher as wm
    atoms = set(atoms)
    out = {"radical_sites": {}, "spin": [], "resonance": []}
    for side in ("reactants", "products"):
        rows = dict(_side_rows(witness, side))
        out["radical_sites"][side] = sorted([a, rows[a][3]] for a in atoms if rows[a][3])
        for part in rm.split_rows(rows):
            if not set(part) & atoms:
                continue
            mol = rm.molecule_from_rows(part)
            verdict = (policy or LEGACY_V1).domain.spin_verdict(mol)
            out["spin"].append({"side": side, "first_atom": min(part), "status": verdict.status,
                                "reason": verdict.reason,
                                "radical_electrons": verdict.radical_electrons})

            def value(a, part=part):
                row = part[a]
                return row[0], row[2], row[4], row[3], dict(row[5])
            try:
                _subs, total = wm._component_subs(sorted(part), value)
                out["resonance"].append({"side": side, "first_atom": min(part), "forms": total})
            except wm.ConjugatedRingError as error:
                out["resonance"].append({"side": side, "first_atom": min(part),
                                         "cutoff": f"conjugated ring > "
                                                   f"{wm.CONJUGATED_RING_LIMIT}: {error}"})
            except wm.ResonanceLimitError:
                out["resonance"].append({"side": side, "first_atom": min(part),
                                         "cutoff": f"> {wm.RESONANCE_FORM_LIMIT} forms"})
    return out


# ---------------------------------------------------------------------------
# record members on provider witnesses
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class MemberCase:
    """One member of a channel on a provider witness: ``witness`` the
    representative mapping (record direction, unextended), ``moves`` every
    mapping of the member, ``wanted`` their family-template label sets."""

    label: str
    family: str
    direction: str
    pathway: tuple
    witness: Any = field(compare=False)
    moves: tuple = field(compare=False)
    wanted: frozenset = field(compare=False)
    support: str = "supported"
    species: tuple = ()
    reactants: tuple = field(default=(), compare=False)
    forms: tuple = field(default=(), compare=False)
    kind: str = "kept"


def template_labels(move, direction) -> frozenset:
    if direction == "forward":
        return frozenset((lab, a) for _s, lab, a in move.labels)
    return frozenset(tuple(x) for x in (move.product_labels or ()))


def channel_moves(certifier, family, molecules, direction, resonance="all"):
    """(reactant molecules the moves are on, their forms, channel set)."""
    from rmgpy.kmc import channels as chm
    from rmgpy.kmc import rmg_moves as rm
    if direction == "forward":
        cs = certifier.channel_set(family, molecules)
        forms = tuple(certifier._forms(cs, k) for k in range(len(cs.reactants)))
        return tuple(cs.reactants), forms, cs
    cs = chm.build_channels(family, molecules, "reverse", resonance=resonance)
    forms = tuple(rm.reactant_forms(m, resonance) for m in cs.reactants)
    return tuple(cs.reactants), forms, cs


def _member_moves(cs, channel):
    from rmgpy.kmc import channels as chm
    from rmgpy.kmc import rmg_moves as rm
    by_outcome = defaultdict(list)
    for m in cs.moves:
        if m.status in chm.ALLOWED_STATUSES and m.template is not None and \
                tuple(m.template) == channel.pathway and m.products is not None:
            by_outcome[rm.product_key(m)].append(m)
    return {o: sorted(v, key=rm.Move.sort_key) for o, v in by_outcome.items()}


def ring_atoms(molecules) -> set:
    out = set()
    for m in molecules:
        out |= {a.id for a in m.atoms if m.is_atom_in_cycle(a)}
    return out


def ring_change_select(family, molecules):
    """``select`` dropping every mapping whose connectivity change (a bond
    formed or broken, ``connectivity_change``) touches any ring atom of the
    reactants.  A generation filter owned by the chemistry policy
    (``ChemistryPolicy.ring_edits == "exclude"``): ``member_cases`` and
    ``probe_slots`` apply it under such a policy and never otherwise.  The
    chemistry it removes is counted (``stats["ring_change_mappings"]``)."""
    ring = ring_atoms(molecules)

    def select(candidate):
        change = connectivity_change(family, _candidate_labels(candidate))
        return not any(a in ring for _kind, pair in change for a in pair)
    return select


def member_cases(certifier, family, builds, direction, *, label="", molecules=None,
                 participants=None, select=None, skip_ring_changes: Optional[bool] = None,
                 stats=None, policy=None):
    """Every member of every channel of ``family`` (``direction``) on the
    provider builds (or given molecules + participants).  With
    ``skip_ring_changes`` mappings whose bond changes touch a ring atom are
    not run through the recipe (``ring_change_select``); ``stats`` counts
    them.  ``skip_ring_changes`` None: the policy's ``ring_edits``
    (``policy``, else ``policy_for`` the witness); an explicit bool keeps its
    meaning."""
    from rmgpy.kmc import certificates as cm
    from rmgpy.kmc import channels as chm
    from rmgpy.kmc import rmg_moves as rm
    if molecules is None:
        molecules, participants = b.witness_reactants(builds)
    if skip_ring_changes is None:
        skip_ring_changes = _builds_policy(builds, policy, participants).ring_edits == \
            "exclude"
    roots = [cm.aromatic_root(m) for m in molecules] if direction == "forward" else \
        list(molecules)
    gate = ring_change_select(family, roots) if skip_ring_changes else None
    cs = chm.build_channels(family, roots, direction, resonance=certifier.resonance,
                            select=gate)
    reactants = tuple(cs.reactants)
    forms = tuple(rm.reactant_forms(m, certifier.resonance) for m in reactants)
    if stats is not None:
        stats["ring_change_mappings"] = stats.get("ring_change_mappings", 0) + sum(
            1 for m in cs.moves if m.status == "excluded")
        stats["mappings"] = stats.get("mappings", 0) + len(cs.moves)
        stats["excluded_outcomes"] = stats.get("excluded_outcomes", 0) + len(cs.excluded)
    out = []
    for ch in cs.channels:
        moves_of = _member_moves(cs, ch)
        for k, outcome in enumerate(ch.members):
            moves = tuple(moves_of.get(outcome, ()))
            if not moves:
                continue
            if select is not None and not select(ch, moves):
                continue
            witness = b.witness_from_move(moves[0], reactants, participants, forms=list(forms))
            species = tuple(ch.reactants) + tuple(ch.products)
            out.append(MemberCase(f"{label}|{family.label}|{direction}|{';'.join(ch.pathway)}|"
                                  f"{','.join(ch.products)}|m{k}", family.label, direction,
                                  tuple(ch.pathway), witness, moves,
                                  frozenset(template_labels(m, direction) for m in moves),
                                  ch.support, species, reactants, forms, ch.kind))
    return out


def species_domain_refusal(witness, policy=None) -> Optional[str]:
    """The policy domain's ``check`` (``policy``, else ``policy_for`` the
    witness) of every species on both sides (None: all accepted)."""
    from rmgpy.kmc import rmg_moves as rm
    domain = (policy or policy_for(witness)).domain
    for side in ("reactants", "products"):
        for part in rm.split_rows(dict(_side_rows(witness, side))):
            verdict = domain.check(rm.molecule_from_rows(part))
            if not verdict.accepted:
                return f"refused:domain: {side} species at atom {min(part)}: {verdict.reason}"
    return None


def mapping_problems(case, window_atoms, fp_witness=None) -> list:
    """Every mapping of the member: reacting atoms (labels, mapped atoms) and
    its edits on BOTH sides inside the window."""
    problems = []
    atoms = set(window_atoms)
    for move in case.moves:
        mw = b.witness_from_move(move, case.reactants, case.witness.participants,
                                 forms=list(case.forms))
        edited = set(b.edit_support(mw))
        labelled = {a for _s, _l, a in move.labels} | {a for _s, _i, _j, a in move.mapped}
        if move.product_labels:
            labelled |= {a for _l, a in move.product_labels}
        outside = sorted((edited | labelled) - atoms)
        if outside:
            problems.append(f"mapping {move.order}/{move.forms}: atoms {outside[:6]} outside "
                            "the window")
    return problems


def run_context(certifier, family, case: MemberCase, assignment, *, bound=None,
                providers=None, start_depths=None, policy=None, **loop) -> ContextResult:
    """A physical context: ends built by the provider on the ORIGINAL witness
    (the end sits where the witness declares it), the event rematched by
    ``rmg_moves`` (same labels, same mapped outcome), the other artificial
    ports pre-extended to ``start_depths`` (the base fixed point's), then
    recertified and closed to its own fixed point."""
    from rmgpy.kmc import rmg_moves as rm
    try:
        built = reconstruct(case.witness, assignment, providers)
    except b.BoundaryError as error:
        return ContextResult(assignment, "refused:end", str(error))
    policy = policy or policy_for(case.witness, providers)
    why = species_domain_refusal(built, policy)
    if why:
        return ContextResult(assignment, "refused:domain", why, built)
    rows = dict(built.reactant_rows())
    molecules = [rm.molecule_from_rows({a: rows[a] for a in p.atom_ids})
                 for p in built.participants]
    rep = case.moves[0]
    try:
        if case.direction == "forward":
            scoped = scoped_certifier(certifier, family, case.wanted)
            reactants, forms, cs = channel_moves(scoped, family, molecules, "forward")
        else:
            from rmgpy.kmc import channels as chm
            cs = chm.build_channels(family, molecules, "reverse",
                                    select=lambda c: c.labels == rep.labels)
            reactants = tuple(cs.reactants)
            forms = tuple(rm.reactant_forms(m) for m in reactants)
    except Exception as error:   # reported, never hidden
        return ContextResult(assignment, "refused:rematch", f"{type(error).__name__}: {error}",
                             built)
    target = rm.product_key(dict(built.product_rows()))
    found = None
    for ch in cs.channels:
        for outcome, moves in _member_moves(cs, ch).items():
            if outcome == target and any(m.labels == rep.labels for m in moves):
                found = (ch, moves)
                break
        if found:
            break
    if found is None:
        same = [m for m in cs.moves if m.labels == rep.labels]
        return ContextResult(assignment, "refused:rematch",
                             f"no move with the labels and the transported outcome "
                             f"({len(same)} moves with the labels: "
                             f"{sorted({m.status for m in same})})", built)
    _ch, moves = found
    move = next(m for m in moves if m.labels == rep.labels)
    witness = b.witness_from_move(move, reactants, built.participants, forms=list(forms))
    wanted = frozenset(template_labels(m, case.direction) for m in moves)
    for port, depth in (start_depths or ()):
        if witness.port(port).artificial:
            for _k in range(depth):
                try:
                    witness, _record = b.extend(witness, port, providers)
                except b.BoundaryError:
                    break
    observer = MemberObserver(certifier, family, case.direction, wanted, providers=providers,
                              policy=policy)
    fp = fixed_point(witness, observer, providers, **loop)
    if fp.reached and bound is not None:
        check_span(fp.witness, fp.window.atoms, bound, providers)
    return ContextResult(assignment, fp.status, fp.reason, fp.witness, fp.window,
                         fp.observation, fp)


def near_ports(witness, atoms, bound) -> tuple:
    """Artificial ports whose port atom lies within ``bound.reach + 1`` heavy
    bonds of ``atoms`` (all ports without a bound)."""
    if bound is None:
        return tuple(b.vector_ports(witness))
    rows = witness.reactant_rows()
    heavy = {a for a, row in rows.items() if row[0] != "H"}
    start = {a for a in atoms if a in heavy}
    start |= {n for a in atoms if a in rows and a not in heavy for n, _o in rows[a][5]
              if n in heavy}
    dist = {a: 0 for a in start}
    frontier = list(dist)
    depth = 0
    while frontier and depth < bound.reach + 1:
        depth += 1
        nxt = []
        for a in frontier:
            for n, _o in rows[a][5]:
                if n in heavy and n not in dist:
                    dist[n] = depth
                    nxt.append(n)
        frontier = nxt
    return tuple(p for p in b.vector_ports(witness) if witness.port(p).atom_id in dist)


def reacting_atoms(case) -> set:
    """The member's reacting atoms on its witness: the edit support and every
    labelled or mapped atom of every mapping."""
    out = set(b.edit_support(case.witness))
    for move in case.moves:
        out |= {a for _s, _l, a in move.labels}
        out |= {a for _s, _i, _j, a in move.mapped}
    return out


def member_near_ports(case, bound) -> tuple:
    """Ports that can lie within the span bound of the member's reacting atoms
    (``near_ports`` around ``reacting_atoms``): the window of any context
    holds atoms within ``bound.reach`` bonds of a reacting atom only, so a
    port farther than ``reach + 1`` heavy bonds (its end atoms, caps and the
    alkene partner included) is outside every window."""
    return near_ports(case.witness, reacting_atoms(case), bound)


# ---------------------------------------------------------------------------
# reacting-port covers (policy a)
# ---------------------------------------------------------------------------
INTERIOR = "interior"
COVER_VERDICTS = ("covered", "out_of_domain", "no_event", "gap")


def reacting_ports(witness) -> tuple:
    """Artificial ports whose port atom is a reacting atom (edit support)."""
    support = b.edit_support(witness)
    return tuple(p for p in b.vector_ports(witness) if witness.port(p).atom_id in support)


def builds_of(specs):
    """Provider builds of witness specs (``build_spec``)."""
    specs = list(specs)
    for spec in specs:
        witness_spec(spec)
    return [build_spec(spec, f"p{k}") for k, spec in enumerate(specs)]


def _longer_declared(spec: DeclaredSpec, head: bool, k: int) -> DeclaredSpec:
    """A declared chain ``k`` periods longer at one end: a head move shifts
    every edit's unit (and a bond edit's other unit) by ``k`` and prepends
    ``k`` copies of the first orientation; a tail move appends ``k`` copies of
    the last.  An end cut at a phase cannot move (``ExtensionRefused``)."""
    from dataclasses import replace
    chain = spec.chain
    phase = chain.head_phase if head else chain.tail_phase
    if phase:
        raise b.ExtensionRefused(f"the {'head' if head else 'tail'} of {spec.provider} is cut "
                                 f"at phase {phase}: no whole period to move")
    orientation = chain.orientation
    edits = chain.edits
    if head:
        if orientation:
            orientation = (orientation[0],) * k + orientation
        edits = tuple(replace(e, unit=e.unit + k,
                              other=None if e.other is None else (e.other[0] + k, e.other[1]))
                      for e in edits)
    elif orientation:
        orientation = orientation + (orientation[-1],) * k
    return DeclaredSpec(spec.provider, replace(chain, units=chain.units + k,
                                               orientation=orientation, edits=edits))


def longer_specs(specs, ports, k: int = 1):
    """``specs`` with every port in ``ports`` (participant, port id) moved
    ``k`` repeats outward (the chain is longer there; the port atom becomes
    interior), and the unit shift of each participant's old coordinates."""
    from dataclasses import replace
    specs, shift = list(specs), defaultdict(int)
    for spec in specs:
        witness_spec(spec)
    for index, port_id in ports:
        spec = specs[index]
        head = port_id == "head"
        if isinstance(spec, DeclaredSpec):
            specs[index] = _longer_declared(spec, head, k)
            if head:
                shift[index] += k
            continue
        if isinstance(spec, tuple) and spec[0] == "alkane":
            specs[index] = ("alkane", spec[1] + k) + tuple(spec[2:])
        else:
            moved = {"units": spec.units + k}
            if head:
                for name in ("radical_units", "secondary_radical_units", "pendant_loss_units"):
                    moved[name] = tuple(u + k for u in getattr(spec, name))
                if spec.ring_radical_unit is not None:
                    moved["ring_radical_unit"] = spec.ring_radical_unit + k
            specs[index] = replace(spec, **moved)
        if head:
            shift[index] += k
    return specs, dict(shift)


def event_specs(specs, assignment):
    """The witness specs of one event a reacting-port member stands for:
    ``assignment`` gives every reacting port either ``INTERIOR`` (the chain one
    repeat longer there: ``longer_specs``) or a declared end chemistry (that
    end built physical there by the provider).  (specs, unit shift); the specs
    are None when a participant's spec cannot carry the end (the alkane
    carries them as ("alkane", n, head end, tail end))."""
    from dataclasses import replace
    longer, shift = longer_specs(specs, [p for p, c in assignment if c == INTERIOR])
    for (index, port_id), choice in assignment:
        if choice == INTERIOR:
            continue
        spec = longer[index]
        if isinstance(spec, DeclaredSpec):
            # a declared chain sets its head / tail end kind; its edits are kept
            if choice not in dc.END_KINDS or choice == "artificial":
                return None, shift
            longer[index] = DeclaredSpec(spec.provider, replace(
                spec.chain, **{"head" if port_id == "head" else "tail": choice}))
            continue
        if isinstance(spec, tuple) and spec[0] == "alkane":
            ends = list(spec[2:4]) + [None] * (4 - len(spec))
            ends[0 if port_id == "head" else 1] = choice
            longer[index] = ("alkane", spec[1], *ends)
            continue
        if not isinstance(spec, b.PsChainSpec):
            return None, shift
        longer[index] = replace(spec, **{"head_end" if port_id == "head" else "tail_end":
                                         choice})
    return longer, shift


def event_signatures(case, shift=None) -> set:
    """Per mapping: {(label, coordinate)} with each labelled atom named by its
    repeat coordinate (participant, unit + shift, position) and a hydrogen by
    its heavy neighbour's coordinate.  The rate rule (pathway) is not part of
    the event: the same mapped event away from the cut point may descend to
    another node (e.g. a methyl H becomes a methylene H)."""
    shift = shift or {}
    coords = {a: (p.index, u + shift.get(p.index, 0), k) for p in case.witness.participants
              for a, u, k in p.metadata.coordinates}
    rows = case.witness.reactant_rows()

    def name(a):
        if a in coords:
            return coords[a]
        heavy = [n for n, _o in rows[a][5] if n in coords]
        return ("H", coords[heavy[0]]) if heavy else ("?", a)
    return {frozenset((lab, name(a)) for _s, lab, a in move.labels) for move in case.moves}


_COVER_CASES = {}


def _cover_cases(certifier, family, specs, direction, policy=None):
    """The member cases of an event's own specs (``member_cases``, the
    generation filter of ``policy`` or of the event witness's policy),
    cached per family, direction, specs and policy name."""
    builds = builds_of(specs)
    policy = _builds_policy(builds, policy)
    key = (family.label, direction, repr(specs), policy.name)
    if key not in _COVER_CASES:
        if len(_COVER_CASES) > 4:
            _COVER_CASES.clear()
        _COVER_CASES[key] = member_cases(certifier, family, builds, direction,
                                         label="cover", policy=policy)
    return _COVER_CASES[key]


def cover_status_verdict(status) -> str:
    """``covered`` (the event's own witness reached a fixed point),
    ``out_of_domain`` (``domain.check`` refuses a species of the event, reason
    stated: ``refused:domain``), ``no_event`` (RMG's family generates no
    mapping on the event's atoms of its own witness, under any labels:
    ``not_generated``; the species has no such reaction, so no rate) or
    ``gap`` (anything else: certificate failures, ceiling hits, incomplete
    descents, port or end refusals, no witness, a rematch that finds the
    atoms under other labels)."""
    status = status or ""
    if status == "fixed_point":
        return "covered"
    if status == "refused:domain":
        return "out_of_domain"
    if status == "not_generated":
        return "no_event"
    return "gap"


def _event_verdict(statuses) -> str:
    verdicts = [cover_status_verdict(s) for s in statuses] or ["gap"]
    if "gap" in verdicts:
        return "gap"
    for verdict in ("covered", "out_of_domain"):
        if verdict in verdicts:
            return verdict
    return "no_event"


def unmatched_event(wanted, built, present) -> tuple:
    """(status, reason) of an event whose own witness has no member with its
    labelled atoms (``wanted``: the event's signatures, ``event_signatures``):
    ``unmatched`` (a gap) when some member (``built``: (label, signatures))
    maps the event's atoms under other labels, or a heavy atom of the event
    is not on the witness (``present``: its repeat coordinates);
    ``not_generated`` (``no_event``, accepted) only when RMG generates no
    mapping on those atoms under any labels."""
    atoms = {frozenset(x for _l, x in sig) for sig in wanted}
    relabelled = [label for label, sigs in built
                  if {frozenset(x for _l, x in sig) for sig in sigs} & atoms]
    heavy = {x for group in atoms for x in group if x[0] != "H"}
    if relabelled or not heavy <= set(present):
        return "unmatched", ("no member has the mapped event "
                             f"({len(relabelled)} map its atoms under other labels)")
    return "not_generated", ("the family generates no mapping on the event's atoms on this "
                             "species (no such event)")


def reacting_port_covers(certifier, family, case, *, bound=None, providers=None, specs=None,
                         policy=None, **loop) -> list:
    """One cover row per distinct mapped event the reacting-port member stands
    for: every assignment of each reacting port to ``INTERIOR`` (the chain
    continues there) or to a declared end chemistry (the chain ends there).
    Each event is run on its OWN witness: the specs with the interior ports one
    repeat longer and the end ports built physical (``event_specs``), the
    member rematched by its labelled atoms' repeat coordinates
    (``event_signatures``; every matching member must be covered), or without
    specs / for an end a spec cannot carry, the end reconstructed on the
    member's witness (``run_context``).  ``verdict``: ``covered``,
    ``out_of_domain`` or ``gap`` (``cover_status_verdict``)."""
    policy = policy or policy_for(case.witness, providers)
    ports = list(b.vector_ports(case.witness))
    reacting = list(reacting_ports(case.witness))
    choices = [(INTERIOR,) + tuple(b.end_chemistries(case.witness, p, providers))
               for p in reacting]
    covers = []
    for combo in product(*choices):
        assignment = tuple(zip(reacting, combo))
        interior = [p for p, c in assignment if c == INTERIOR]
        ends = [c for _p, c in assignment if c != INTERIOR]
        row = {"event": ", ".join(f"{c}@{p[0]}:{p[1]}" for p, c in assignment),
               "kind": "interior" if not ends else "physical" if not interior else "mixed",
               "port": [list(p) for p in reacting] if len(reacting) > 1 else list(reacting[0]),
               "end": None if not ends else ends[0] if len(ends) == 1 else ends,
               "assignment": [[list(p), c] for p, c in assignment]}
        especs, shift = (None, {}) if specs is None else event_specs(specs, assignment)
        built = None
        if especs is not None:
            row.update(route="spec", specs=[repr(x) for x in especs])
            try:
                builds = builds_of(especs)
            except b.BoundaryError as error:
                row.update(status="refused:end", reason=str(error)[:200])
            else:
                refusals = [(k, policy.domain.check(x.molecule)) for k, x in enumerate(builds)]
                refusals = [(k, v) for k, v in refusals if not v.accepted]
                if refusals:
                    # the event's own reactants are outside the declared domain
                    row.update(status="refused:domain", reason="; ".join(
                        f"reactants species {k}: {v.reason}" for k, v in refusals)[:300])
                else:
                    built = _cover_cases(certifier, family, especs, case.direction, policy)
        if built is not None:
            wanted = event_signatures(case, shift)
            matches = [c for c in built if event_signatures(c) & wanted]
            if not matches:
                present = {(p.index, u, k) for p in built[0].witness.participants
                           for _a, u, k in p.metadata.coordinates} if built else \
                    {(i, u, k) for i, x in enumerate(builds) for _a, u, k in
                     x.metadata.coordinates}
                status, reason = unmatched_event(
                    wanted, [(c.label, event_signatures(c)) for c in built], present)
                row.update(status=status, reason=reason)
            else:
                members = []
                for c in matches:
                    # an event with a species outside the declared domain needs no
                    # rate (proved by domain.check, reason stated), whatever else
                    # refuses it first
                    why = species_domain_refusal(c.witness, policy)
                    res = {"status": "refused:domain", "reason": why} if why else \
                        run_member(certifier, family, c, bound=bound, contexts=False,
                                   providers=providers, covers=False, policy=policy, **loop)
                    members.append({"label": c.label[-120:], "pathway": list(c.pathway),
                                    "status": res["status"],
                                    "reason": str(res.get("reason", ""))[:200],
                                    "window_atoms": res.get("window_atoms")})
                row["members"] = members
                row["status"] = ",".join(sorted({m["status"] for m in members}))
                row["reason"] = "; ".join(sorted({m["reason"] for m in members
                                                  if m["reason"]}))[:300]
                row["window_atoms"] = members[0]["window_atoms"]
        elif "status" not in row and not interior:
            r = run_context(certifier, family, case, tuple((p, dict(assignment).get(p))
                                                           for p in ports),
                            bound=bound, providers=providers, policy=policy, **loop)
            row.update(route="reconstruct", status=r.status, reason=r.reason[:200],
                       members=[{"status": r.status, "reason": r.reason[:200]}],
                       window_atoms=len(r.window.atoms) if r.status == "fixed_point" else None)
        elif "status" not in row:
            row.update(status="no_witness", reason="the interior event needs the witness specs")
        statuses = [m["status"] for m in row.get("members", ())] or [row["status"]]
        row["verdict"] = _event_verdict(statuses)
        covers.append(row)
    return covers


def _event_key(cover) -> tuple:
    if cover.get("event") is not None:
        return (cover["event"],)
    return (cover.get("kind"), repr(cover.get("port")), repr(cover.get("end")))


def cover_verdict(covers) -> bool:
    """A reacting-port member is covered iff EVERY distinct mapped event it
    stands for (the interior event and each physical-end event, ``_event_key``)
    is covered by its own witness reaching a fixed point or shown outside the
    declared domain (``cover_status_verdict``); a cover of one event never
    stands in for another, and no event, or any gap, is uncovered."""
    events = defaultdict(list)
    for c in covers:
        events[_event_key(c)].append(c)
    if not events:
        return False
    for rows in events.values():
        for c in rows:
            statuses = [m.get("status") for m in c.get("members", ())] or [c.get("status")]
            if _event_verdict(statuses) == "gap":
                return False
    return True


def cover_gaps(covers) -> list:
    """The events of a cover list that are gaps (``cover_verdict``)."""
    out = []
    for c in covers:
        statuses = [m.get("status") for m in c.get("members", ())] or [c.get("status")]
        if _event_verdict(statuses) == "gap":
            out.append(f"{c.get('event') or (c.get('kind'), c.get('port'), c.get('end'))}: "
                       f"{','.join(sorted({str(s) for s in statuses}))}")
    return out


def reacting_port_gaps(rows) -> list:
    """Every reacting-port member with an uncovered event (``cover_verdict`` on
    its ``covers``; a member without a cover search is a gap)."""
    out = []
    for r in rows:
        if r.get("status") != "refused:reacting_port":
            continue
        covers = r.get("covers") or []
        if not cover_verdict(covers):
            out.append(f"{r.get('label', '?')[-90:]}: reacting port with uncovered events "
                       f"{cover_gaps(covers) or ['no cover search']}")
    return out


def run_member(certifier, family, case: MemberCase, *, bound=None, contexts: bool = True,
               providers=None, specs=None, covers: bool = True, keep=None, contexts_out=None,
               policy=None, **loop) -> dict:
    """Items 1-5 and 7 for one member: refusal by status, or the fixed point
    (window, span against ``bound``, per-mapping edit checks, diagnostics,
    recertification flips) and the context list with merges.  A
    reacting-port refusal carries its covering witnesses (``covers``,
    ``covered``; the interior cover needs the witness ``specs``).
    ``contexts_out`` (additive, I-138): a list that receives the context
    results themselves (``ContextResult`` per assignment, first the
    all-artificial base, and the merge pairs) as one ``(results, merges)``
    entry; the returned rows are unchanged.  ``policy``: the chemistry
    (default ``policy_for`` the member's witness), passed to every
    observation, context and cover."""
    import time
    start = time.perf_counter()
    out = {"label": case.label, "family": case.family, "direction": case.direction,
           "pathway": list(case.pathway), "mappings": len(case.moves),
           "slot_orders": sorted({m.order for m in case.moves}),
           "resonance_variants": sorted({m.forms for m in case.moves})}
    if case.support != "supported":
        out.update(status=f"refused:{case.support}", reason="channel support")
        return out
    policy = policy or policy_for(case.witness, providers)
    why = species_domain_refusal(case.witness, policy)
    if why:
        out.update(status="refused:domain", reason=why)
        return out
    observer = MemberObserver(certifier, family, case.direction, case.wanted,
                              providers=providers, policy=policy)
    fp = fixed_point(case.witness, observer, providers, **loop)
    if keep is not None:       # the class closings keep the fixed point (``close_classes``)
        keep.append(fp)
    out.update(status=fp.status, reason=fp.reason, depths=[[list(p), d] for p, d in fp.depths],
               retries=fp.retries, observations=observer.calls,
               flips=sorted({(n, o, w) for _k, n, o, w in fp.flips}),
               check_rounds=len(fp.checks))
    if not fp.reached:
        if fp.status == "refused:reacting_port" and covers:
            out["covers"] = reacting_port_covers(certifier, family, case, bound=bound,
                                                 providers=providers, specs=specs,
                                                 policy=policy, **loop)
            out["covered"] = cover_verdict(out["covers"])
        out["seconds"] = round(time.perf_counter() - start, 1)
        return out
    window = fp.window
    out["window_atoms"] = len(window.atoms)
    out["window_heavy"] = sum(1 for a in window.atoms if fp.witness.reactant_rows()[a][0] != "H")
    out["boundary_bond_problems"] = boundary_bond_problems(fp.witness, window)
    out["mapping_problems"] = mapping_problems(case, window.atoms)
    span = window_span(fp.witness, window.atoms, providers)
    out["span"] = {str(k): v for k, v in span.items()}
    if bound is not None:
        check_span(fp.witness, window.atoms, bound, providers)
        out["span_bound_held"] = True
    out["dispositions"] = sorted({(str(k[0]), d) for k, d, _a in window.support.dispositions})
    out["diagnostics"] = diagnostics(fp.witness, window.atoms, policy)
    if contexts:
        _contexts(out, certifier, family, case, fp, bound, providers, loop,
                  contexts_out=contexts_out, policy=policy)
    out["seconds"] = round(time.perf_counter() - start, 1)
    return out


def _contexts(out, certifier, family, case, fp, bound, providers, loop, contexts_out=None,
              policy=None):
    """Every context of ``context_assignments`` over the member's near ports
    (``member_near_ports``), run from the base fixed point's depths; far
    ports are reduced to their radical ends, and the reduction is checked:
    each far non-radical end single must merge with the all-artificial
    context, otherwise that port becomes near and the product is rebuilt
    (``far_port_promotions``)."""
    ports = list(b.vector_ports(case.witness))
    base = ContextResult(tuple((p, None) for p in ports), fp.status, fp.reason, fp.witness,
                         fp.window, fp.observation, fp)
    near = set(member_near_ports(case, bound))
    done = {base.assignment: base}
    promoted = []
    while True:
        assignments = context_assignments(case.witness, near, providers)
        for assignment in assignments[1:]:
            if assignment not in done:
                done[assignment] = run_context(certifier, family, case, assignment, bound=bound,
                                               providers=providers, start_depths=fp.depths,
                                               policy=policy, **loop)
        results = [done[a] for a in assignments]
        far = [p for p in ports if p not in near]
        verdicts = far_port_verdicts(results, far, providers, witness=case.witness)
        more = far_port_promotions(verdicts)
        if not more:
            break
        promoted.extend(more)
        near |= set(more)
    out["context_policy"] = {
        "near_ports": [list(p) for p in ports if p in near],
        "far_ports": [list(p) for p in ports if p not in near],
        "promoted": [list(p) for p in promoted],
        "far_checks": [[[list(p) for p in _verdict_ports(k)],
                        [k[1]] if isinstance(k[1], str) else [e for _p, e in k], ok]
                       for k, ok in sorted(verdicts.items(), key=repr)]}
    merges = merge_contexts(results)
    merged = {j for _i, j in merges}
    rows = []
    for k, r in enumerate(results):
        span = window_span(r.witness, r.window.atoms, providers) \
            if r.window and r.status == "fixed_point" \
            else None
        rows.append({"assignment": [[list(p), e] for p, e in r.assignment if e],
                     "status": r.status, "reason": r.reason[:200],
                     "window_atoms": len(r.window.atoms) if span else None,
                     "both_ends": any(s["both_ends"] for s in span.values()) if span else None,
                     "physical_ends_inside": {str(i): s["physical_ends_inside"]
                                              for i, s in span.items()} if span else None,
                     "merged_into": next((i for i, j in merges if j == k), None)})
    out["contexts"] = rows
    out["context_classes"] = len([k for k, r in enumerate(results)
                                  if r.status == "fixed_point" and k not in merged])
    out["merges"] = [list(x) for x in merges]
    if contexts_out is not None:
        contexts_out.append((tuple(results), tuple(merges)))
