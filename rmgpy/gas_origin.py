#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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
Gas-origin provenance for reaction generation.

Rule: a generated molecule whose heavy atoms ALL originate from genuine gas
reactant molecules is gas -- never a Polymer, never a pool daughter, never
folded into a pool by any structural path. Everything else keeps the
structural behaviour.

Mechanics (single generation call, ``KineticsDatabase.generate_reactions_from_families``):

1. Each supplied reactant is classified with :func:`is_genuine_gas_reactant`.
2. Heavy atoms of COPIES of the genuine-gas reactants' molecules are stamped
   ``atom.props[GAS_ORIGIN_KEY] = k`` (k = the reactant's occurrence index in
   the supplied tuple; presence is tested with ``in``, never truthiness,
   because k may be 0). Originals are never touched.
3. After the families run, every reaction gets explicit metadata
   (:func:`assign_gas_origin_metadata`): which side holds the generated
   participants, one gas verdict per generated participant, and one source
   index per supplied participant.
4. The degeneracy collapse refuses merges whose (structure, verdict)
   multisets differ (:func:`gas_origin_multisets_differ`).
5. All stamps are stripped (:func:`strip_gas_origin`) before the reactions
   leave the generation call.

``CoreEdgeReactionModel.make_new_reaction`` consumes the metadata: supplied
participants are restored by source index, gas-verdict participants are
resolved with a lookup that never returns a condensed species.
"""

from collections import Counter

#: Atom ``props`` key carrying the occurrence index of the genuine gas
#: reactant an atom came from. Transient: present only inside
#: ``generate_reactions_from_families``.
GAS_ORIGIN_KEY = "gas_origin"

#: ``polymer_refused_detail`` of a degeneracy-collapse representative whose
#: merged rows disagreed on which participants are gas.
ORIGIN_AMBIGUOUS_COLLAPSE = "origin-ambiguous-collapse"

#: Census keys, one place so every producer and consumer agrees.
CENSUS_COLLAPSE_REFUSALS = "collapse_refusals"
CENSUS_IDENTITY_DROPS = "identity_filter_drops"
CENSUS_IDENTITY_DROPS_GAS_MIXED = "identity_filter_drops_gas_mixed"
CENSUS_GAS_FOLD_BYPASSED = "gas_fold_bypassed"
CENSUS_GAS_LOOKUP_UNPOOLED = "gas_lookup_unpooled"
CENSUS_SUPPLIED_RELOOKUP_DIFFERS = "supplied_relookup_differs"
CENSUS_ROWS_CHANGED_CLASS = "rows_changed_class"

CENSUS_KEYS = (
    CENSUS_COLLAPSE_REFUSALS,
    CENSUS_IDENTITY_DROPS,
    CENSUS_IDENTITY_DROPS_GAS_MIXED,
    CENSUS_GAS_FOLD_BYPASSED,
    CENSUS_GAS_LOOKUP_UNPOOLED,
    CENSUS_SUPPLIED_RELOOKUP_DIFFERS,
    CENSUS_ROWS_CHANGED_CLASS,
)


def new_census():
    """Return an empty gas-origin census (a :class:`collections.Counter`
    holding every key at zero, so reports never omit a key)."""
    census = Counter()
    for key in CENSUS_KEYS:
        census[key] = 0
    return census


def merge_census(target, source):
    """Add the counts of ``source`` into ``target`` (both dict-like)."""
    if target is None or not source:
        return
    for key, value in source.items():
        target[key] = target.get(key, 0) + value


def format_census(census):
    """One-line ``key=value`` rendering in :data:`CENSUS_KEYS` order."""
    return ", ".join("{0}={1}".format(k, int(census.get(k, 0))) for k in CENSUS_KEYS)


################################################################################
# 1. Genuine gas
################################################################################

def is_genuine_gas_reactant(species, condensed=False):
    """Return True if ``species`` is a genuine gas reactant.

    Not genuine gas: a :class:`~rmgpy.polymer.Polymer`; a moment dummy; an
    explicit-DP oligomer (``props['explicit_dp_origin']``, set in
    ``rmg/input.py``); or anything the caller's condensed-membership
    predicate classifies condensed (``condensed=True``). The latter is the
    solver's own core phase mask (``PolymerPhase.get_gas_mask``: configured
    condensed proxies, explicit-DP members, condensed H-loss daughters,
    everything booked into a pool), evaluated by the caller over the full
    species list and passed in here -- this function never re-derives it and
    applies no size rule.
    """
    if condensed:
        return False
    from rmgpy.polymer import Polymer
    if isinstance(species, Polymer) or getattr(species, "is_polymer", False):
        return False
    if getattr(species, "is_moment_dummy", False):
        return False
    props = getattr(species, "props", None)
    if isinstance(props, dict) and "explicit_dp_origin" in props:
        return False
    return True


################################################################################
# 2. Stamp
################################################################################

def stamp_gas_origin_copy(molecule, k):
    """Return a deep copy of ``molecule`` whose heavy atoms carry
    ``props[GAS_ORIGIN_KEY] = k``. Atom ids are preserved by the copy, so the
    degeneracy bookkeeping (identical-reaction detection) is unaffected."""
    stamped = molecule.copy(deep=True)
    for atom in stamped.atoms:
        if not atom.is_hydrogen():
            atom.props[GAS_ORIGIN_KEY] = k
    return stamped


def _participant_molecules(participant):
    """Every molecule graph held by a participant (a Molecule, or every
    resonance structure of a Species/Polymer)."""
    mols = getattr(participant, "molecule", None)
    if mols is not None and not hasattr(participant, "atoms"):
        return [m for m in mols if m is not None]
    return [participant]


def _first_molecule(participant):
    mols = _participant_molecules(participant)
    return mols[0] if mols else None


def molecule_is_gas_origin(molecule):
    """True iff ``molecule`` has at least one heavy atom and EVERY heavy atom
    carries a gas-origin stamp (key presence: the stamp value may be 0).
    Heavy-atom-free molecules carry no provenance and return False (they keep
    the structural behaviour)."""
    if molecule is None:
        return False
    heavy = [a for a in molecule.atoms if not a.is_hydrogen()]
    if not heavy:
        return False
    return all(GAS_ORIGIN_KEY in a.props for a in heavy)


################################################################################
# 3. Verdict metadata
################################################################################

def generated_side(rxn):
    """Name of the side holding the GENERATED participants of a
    family-generated reaction: ``'products'`` for a forward-written row,
    ``'reactants'`` for a reverse-written one (``family._create_reaction``
    puts the generated structures in ``reactants`` when ``is_forward`` is
    False)."""
    return "products" if rxn.is_forward else "reactants"


def _other_side(side):
    return "reactants" if side == "products" else "products"


def assign_gas_origin_metadata(rxn, combo):
    """Stamp explicit gas-origin metadata onto a freshly generated reaction
    whose participants are still the raw Molecule objects of ``combo`` (the
    supplied side) and the generated product structures.

    * ``rxn.gas_origin_side`` -- the generated side name.
    * ``rxn.gas_origin_verdicts`` -- one bool per participant on that side.
    * ``rxn.gas_origin_sources`` -- one supplied-reactant index (or None when
      the participant is not one of the combo objects) per participant on
      the other side.
    """
    side = generated_side(rxn)
    generated = getattr(rxn, side)
    supplied = getattr(rxn, _other_side(side))
    rxn.gas_origin_side = side
    rxn.gas_origin_verdicts = tuple(
        molecule_is_gas_origin(_first_molecule(p)) for p in generated)
    sources = []
    for participant in supplied:
        src = None
        for k, mol in enumerate(combo):
            if participant is mol:
                src = k
                break
        sources.append(src)
    rxn.gas_origin_sources = tuple(sources)


def restore_supplied_originals(rxn, originals):
    """Swap the stamped COPIES on a generated row's supplied side (and in its
    ``pairs``) back to the original, never-stamped molecule objects, so the
    row's object graph is exactly what generation produced before
    provenance existed. ``originals`` maps ``id(copy) -> original``."""
    if not originals:
        return
    side = _other_side(generated_side(rxn))
    current = getattr(rxn, side)
    setattr(rxn, side, [originals.get(id(p), p) for p in current])
    if rxn.pairs:
        rxn.pairs = [type(pair)((originals.get(id(pair[0]), pair[0]), originals.get(id(pair[1]), pair[1])))
                     for pair in rxn.pairs]


def has_gas_origin_metadata(rxn):
    """True if ``rxn`` carries complete, still-aligned gas-origin metadata."""
    side = getattr(rxn, "gas_origin_side", None)
    verdicts = getattr(rxn, "gas_origin_verdicts", None)
    sources = getattr(rxn, "gas_origin_sources", None)
    if side not in ("reactants", "products") or verdicts is None or sources is None:
        return False
    return (len(verdicts) == len(getattr(rxn, side))
            and len(sources) == len(getattr(rxn, _other_side(side))))


def copy_gas_origin_metadata(source, target):
    """Copy the metadata attributes from ``source`` onto ``target``."""
    target.gas_origin_side = getattr(source, "gas_origin_side", None)
    target.gas_origin_verdicts = getattr(source, "gas_origin_verdicts", None)
    target.gas_origin_sources = getattr(source, "gas_origin_sources", None)


def flip_gas_origin_metadata(rxn):
    """Keep the metadata honest after an in-place reactant/product swap: the
    verdict and source tuples move with their lists, only the side name
    changes."""
    side = getattr(rxn, "gas_origin_side", None)
    if side in ("reactants", "products"):
        rxn.gas_origin_side = _other_side(side)


def capture_gas_origin_alignment(rxn):
    """Snapshot (participant, value) pairs so the metadata can be realigned
    after the participant lists are reordered in place (e.g. sorted).
    Returns None when there is no metadata."""
    if not has_gas_origin_metadata(rxn):
        return None
    side = rxn.gas_origin_side
    return (side,
            list(zip(getattr(rxn, side), rxn.gas_origin_verdicts)),
            list(zip(getattr(rxn, _other_side(side)), rxn.gas_origin_sources)))


def restore_gas_origin_alignment(rxn, snapshot):
    """Realign metadata captured by :func:`capture_gas_origin_alignment` to
    the (reordered) participant lists, matching participants by identity."""
    if snapshot is None:
        return
    side, gen_pairs, sup_pairs = snapshot

    def _realign(current, pairs):
        remaining = list(pairs)
        out = []
        for participant in current:
            for j, (obj, value) in enumerate(remaining):
                if obj is participant:
                    out.append(value)
                    del remaining[j]
                    break
            else:
                return None
        return tuple(out)

    verdicts = _realign(getattr(rxn, side), gen_pairs)
    sources = _realign(getattr(rxn, _other_side(side)), sup_pairs)
    if verdicts is None or sources is None:
        return
    rxn.gas_origin_side = side
    rxn.gas_origin_verdicts = verdicts
    rxn.gas_origin_sources = sources


################################################################################
# Strip
################################################################################

def strip_gas_origin(rxn):
    """Remove every gas-origin stamp reachable from ``rxn``: both sides, every
    resonance structure of every participant, and the same for the reverse
    reaction retained by own-reverse families (``rxn.reverse``)."""
    for reaction in (rxn, getattr(rxn, "reverse", None)):
        if reaction is None:
            continue
        for participant in list(reaction.reactants or []) + list(reaction.products or []):
            for mol in _participant_molecules(participant):
                for atom in mol.atoms:
                    atom.props.pop(GAS_ORIGIN_KEY, None)


def any_gas_origin_stamp(rxn):
    """True if any gas-origin stamp is reachable from ``rxn`` (see
    :func:`strip_gas_origin` for the reach)."""
    for reaction in (rxn, getattr(rxn, "reverse", None)):
        if reaction is None:
            continue
        for participant in list(reaction.reactants or []) + list(reaction.products or []):
            for mol in _participant_molecules(participant):
                if any(GAS_ORIGIN_KEY in atom.props for atom in mol.atoms):
                    return True
    return False


################################################################################
# 4. Collapse
################################################################################

def _generated_pairs(rxn):
    if not has_gas_origin_metadata(rxn):
        return None
    return list(zip(getattr(rxn, rxn.gas_origin_side), rxn.gas_origin_verdicts))


def _same_structure(a, b):
    if hasattr(a, "molecule") and not hasattr(a, "atoms"):
        return a.is_isomorphic(b, strict=False)
    if hasattr(b, "molecule") and not hasattr(b, "atoms"):
        return b.is_isomorphic(a, strict=False)
    return a.is_isomorphic(b, strict=False)


def gas_origin_multisets_differ(rxn_a, rxn_b):
    """True iff the multisets of (generated participant structure, gas
    verdict) of two structurally isomorphic generated rows differ. Rows
    without metadata never differ (legacy behaviour)."""
    pairs_a = _generated_pairs(rxn_a)
    pairs_b = _generated_pairs(rxn_b)
    if pairs_a is None or pairs_b is None:
        return False
    if len(pairs_a) != len(pairs_b):
        return True
    remaining = list(pairs_b)
    for part, verdict in pairs_a:
        for j, (other, other_verdict) in enumerate(remaining):
            if other_verdict == verdict and _same_structure(part, other):
                del remaining[j]
                break
        else:
            return True
    return False


def mark_origin_ambiguous_collapse(rxn):
    """Refuse the retained representative of an origin-ambiguous merge (OR
    semantics: an existing refusal and its detail are kept)."""
    rxn.polymer_refused = True
    if getattr(rxn, "polymer_refused_detail", None) is None:
        rxn.polymer_refused_detail = ORIGIN_AMBIGUOUS_COLLAPSE


################################################################################
# Identity filter accounting
################################################################################

def participants_mix_gas_origin(participants):
    """True if the participant list holds at least one fully gas-stamped
    molecule AND at least one heavy-atom-bearing molecule that is not --
    i.e. gas and non-gas provenance meet in the row."""
    has_gas = has_other = False
    for participant in participants:
        mol = _first_molecule(participant)
        if mol is None or not any(not a.is_hydrogen() for a in mol.atoms):
            continue
        if molecule_is_gas_origin(mol):
            has_gas = True
        else:
            has_other = True
    return has_gas and has_other
