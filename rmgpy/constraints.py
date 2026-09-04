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

import logging
import os

from collections import Counter

from rmgpy.species import Species


# ---------------------------------------------------------------------------
# Polymer-constraint warn-once state
# ---------------------------------------------------------------------------

_unbounded_polymer_warned = False

# ---------------------------------------------------------------------------
# Generation census (I-067)
# ---------------------------------------------------------------------------
#
# A size bound set too low truncates the mechanism silently: the run gets fast,
# the log looks clean, and the chemistry that would have mattered was simply
# never generated. Nothing in RMG's output announces it. This census does --
# every structure a constraint tier refuses is counted, bucketed by heavy-atom
# count and by which constraint refused it, and up to _CENSUS_SMILES_CAP distinct
# refused structures are kept by SMILES so the cut set can be named and
# cross-referenced against an unbounded run's edge fluxes.
#
# Only the REFUSAL path is instrumented. A refusal raises ForbiddenStructureException
# at every call site, which dwarfs the cost of a to_smiles(); the accept path --
# the genuinely hot one -- is untouched. The size distribution of what was
# RETAINED is read off the finished model's species lists instead, which is both
# cheaper and exact.

_CENSUS_SMILES_CAP = 20000

# Under RMG_I067_PROFILE the census stops rendering SMILES. The wasted-build
# profile times the refusal path, and to_smiles() on a refused structure is this
# census's own cost, not the generator's -- leaving it in would inflate the very
# fraction the profile exists to measure. Counts are unaffected.
_CENSUS_SMILES_OFF = os.environ.get('RMG_I067_PROFILE') == '1'

_census = {}


def reset_generation_census():
    """Clear the constraint-refusal census. Called from RMG.initialize."""
    _census.clear()
    _census.update({
        'refused': Counter(),                                    # tier -> n
        'by_heavy': {'gas': Counter(), 'polymer': Counter()},    # tier -> heavy -> n
        'by_reason': Counter(),                                  # 'tier/constraint' -> n
        'smiles': {'gas': {}, 'polymer': {}},                    # tier -> smiles -> [n, heavy, C]
        'smiles_capped': False,
    })


reset_generation_census()


def get_generation_census():
    """Return the live constraint-refusal census dict (not a copy)."""
    return _census


def _record_refusal(struct, tier, reason):
    """Record one constraint refusal. Reached only when a structure is rejected."""
    try:
        heavy = struct.get_num_atoms() - struct.get_num_atoms('H')
        carbon = struct.get_num_atoms('C')
    except Exception:
        return
    _census['refused'][tier] += 1
    _census['by_heavy'][tier][heavy] += 1
    _census['by_reason']['{0}/{1}'.format(tier, reason.split(':')[0])] += 1
    if _CENSUS_SMILES_OFF:
        return
    bucket = _census['smiles'][tier]
    if len(bucket) >= _CENSUS_SMILES_CAP:
        _census['smiles_capped'] = True
        return
    try:
        smi = struct.to_smiles()
    except Exception:
        return
    entry = bucket.get(smi)
    if entry is None:
        bucket[smi] = [1, heavy, carbon]
    else:
        entry[0] += 1


def log_generation_census(header='GENERATION CENSUS'):
    """
    Log what the constraint tiers refused: by tier, by heavy-atom count and by
    constraint. Emitted after every enlargement, so the truncation a bound causes
    is visible while it happens rather than only in a post-hoc diff of two runs.
    """
    total = sum(_census['refused'].values())
    if not total:
        logging.info('%s: no structure was refused by any constraint tier.', header)
        return
    logging.info('%s: %d structure-refusals (gas tier %d, polymer tier %d)%s',
                 header, total, _census['refused']['gas'], _census['refused']['polymer'],
                 '; SMILES capture CAPPED at {0} distinct'.format(_CENSUS_SMILES_CAP)
                 if _census['smiles_capped'] else '')
    for tier in ('gas', 'polymer'):
        hist = _census['by_heavy'][tier]
        if not hist:
            continue
        logging.info('    %s tier refusals by heavy-atom count: %s', tier,
                     ', '.join('{0}:{1}'.format(h, hist[h]) for h in sorted(hist)))
        logging.info('    %s tier distinct refused structures: %d', tier,
                     len(_census['smiles'][tier]))
    for key in sorted(_census['by_reason']):
        logging.info('    refused by %s: %d', key, _census['by_reason'][key])

# Default heavy-atom count at or above which a discrete species routes to the polymer
# tier (generatePolymerConstraints) rather than the gas tier. Polymer-independent (a flat
# count), overridable per deck via generatePolymerConstraints(polymerSizeThreshold=...).
DEFAULT_POLYMER_SIZE_THRESHOLD = 15


def reset_polymer_warning():
    """Reset the once-per-run unbounded-polymer warning flag. Called from RMG.initialize."""
    global _unbounded_polymer_warned
    _unbounded_polymer_warned = False


def _warn_unbounded_polymer_once():
    global _unbounded_polymer_warned
    if not _unbounded_polymer_warned:
        # Scope note (I-067): this fires ONLY for real Polymer objects. Generated
        # proxy structures are Molecules/Species, which carry is_polymer_proxy but
        # never is_polymer, and with polymer_constraints=None size routing is off
        # (is_polymer_constraint_member returns False), so they are bounded by the
        # gas tier like everything else. The message used to say "polymer proxy
        # chemistry is unbounded", which was false and sent readers looking for a
        # runaway generation leg that does not exist.
        logging.warning(
            "Polymer pool objects bypass generatedSpeciesConstraints because "
            "generatePolymerConstraints is not set. Generated proxy species are NOT "
            "affected -- with no polymer tier configured they fall through to "
            "generatedSpeciesConstraints. Adding a generatePolymerConstraints block "
            "switches on size routing: every species with heavy-atom count >= "
            "polymerSizeThreshold (default {0}) then leaves the gas tier for the "
            "polymer tier, so a block looser than generatedSpeciesConstraints "
            "ENLARGES the generated space rather than bounding it.".format(
                DEFAULT_POLYMER_SIZE_THRESHOLD)
        )
        _unbounded_polymer_warned = True


def is_polymer_constraint_member(obj, polymer_constraints=None):
    """
    True if `obj` (a Species, Polymer, or Molecule) should be bounded by the polymer tier
    (generatePolymerConstraints) rather than the global generatedSpeciesConstraints.

    Routing is provenance-by-size, not by the leaky is_polymer_proxy flag:
      * real Polymer objects carry the durable `is_polymer` flag -> always polymer tier;
      * any other species is a polymer member iff it is polymer-sized -- heavy-atom count at
        or above `polymerSizeThreshold` (default DEFAULT_POLYMER_SIZE_THRESHOLD).
    The heavy-atom threshold is polymer-independent (a flat count), so the same knob works for
    PE, PP, novolak, etc. Size routing needs a configured polymer tier; when
    `polymer_constraints` is None only real Polymer objects route to polymer. Duck-typed to
    avoid a Polymer import / circular dependency.
    """
    if getattr(obj, 'is_polymer', False):
        return True
    if polymer_constraints is None:
        return False
    threshold = polymer_constraints.get('polymerSizeThreshold', DEFAULT_POLYMER_SIZE_THRESHOLD)
    struct = _normalize(obj)
    heavy = struct.get_num_atoms() - struct.get_num_atoms('H')
    return heavy >= threshold


def _normalize(species):
    """Return the Molecule to constrain: molecule[0] for a Species, else the object itself."""
    if isinstance(species, Species):
        return species.molecule[0]
    return species


def _evaluate_constraints(struct, constraints):
    """Return a reason string if `struct` violates `constraints`, else False. Pure: no get_input."""
    explicitly_allowed_molecules = constraints.get('explicitlyAllowedMolecules', [])
    for molecule in explicitly_allowed_molecules:
        if struct.is_isomorphic(molecule):
            return False

    max_carbon_atoms = constraints.get('maximumCarbonAtoms', -1)
    if max_carbon_atoms != -1:
        if struct.get_num_atoms('C') > max_carbon_atoms:
            return f"Exceeded maximumCarbonAtoms: {struct.get_num_atoms('C')} > {max_carbon_atoms}"

    max_oxygen_atoms = constraints.get('maximumOxygenAtoms', -1)
    if max_oxygen_atoms != -1:
        if struct.get_num_atoms('O') > max_oxygen_atoms:
            return f"Exceeded maximumOxygenAtoms: {struct.get_num_atoms('O')} > {max_oxygen_atoms}"

    max_nitrogen_atoms = constraints.get('maximumNitrogenAtoms', -1)
    if max_nitrogen_atoms != -1:
        if struct.get_num_atoms('N') > max_nitrogen_atoms:
            return f"Exceeded maximumNitrogenAtoms: {struct.get_num_atoms('N')} > {max_nitrogen_atoms}"

    max_silicon_atoms = constraints.get('maximumSiliconAtoms', -1)
    if max_silicon_atoms != -1:
        if struct.get_num_atoms('Si') > max_silicon_atoms:
            return f"Exceeded maximumSiliconAtoms: {struct.get_num_atoms('Si')} > {max_silicon_atoms}"

    max_sulfur_atoms = constraints.get('maximumSulfurAtoms', -1)
    if max_sulfur_atoms != -1:
        if struct.get_num_atoms('S') > max_sulfur_atoms:
            return f"Exceeded maximumSulfurAtoms: {struct.get_num_atoms('S')} > {max_sulfur_atoms}"

    max_heavy_atoms = constraints.get('maximumHeavyAtoms', -1)
    if max_heavy_atoms != -1:
        heavy_atoms = struct.get_num_atoms() - struct.get_num_atoms('H')
        if heavy_atoms > max_heavy_atoms:
            return f"Exceeded maximumHeavyAtoms: {heavy_atoms} > {max_heavy_atoms}"

    max_surface_sites = constraints.get('maximumSurfaceSites', -1)
    if max_surface_sites != -1:
        if struct.get_num_atoms('X') > max_surface_sites:
            return f"Exceeded maximumSurfaceSites: {struct.get_num_atoms('X')} > {max_surface_sites}"

    max_surface_bond_order = constraints.get('maximumSurfaceBondOrder', -1)
    if max_surface_bond_order != -1:
        for site in struct.get_surface_sites():
            if site.get_total_bond_order() > max_surface_bond_order:
                return f"Exceeded maximumSurfaceBondOrder at site: {site.get_total_bond_order()} > {max_surface_bond_order}"

    max_radicals = constraints.get('maximumRadicalElectrons', -1)
    if max_radicals != -1:
        if struct.get_radical_count() > max_radicals:
            return f"Exceeded maximumRadicalElectrons: {struct.get_radical_count()} > {max_radicals}"

    max_carbenes = constraints.get('maximumSingletCarbenes', 1)
    if max_carbenes != -1:
        if struct.get_singlet_carbene_count() > max_carbenes:
            return f"Exceeded maximumSingletCarbenes: {struct.get_singlet_carbene_count()} > {max_carbenes}"

    max_carbene_radicals = constraints.get('maximumCarbeneRadicals', 0)
    if max_carbene_radicals != -1:
        if struct.get_singlet_carbene_count() > 0 and struct.get_radical_count() > max_carbene_radicals:
            return f"Exceeded maximumCarbeneRadicals: {struct.get_radical_count()} > {max_carbene_radicals}"

    return False


def pass_cutting_threshold(species):
    """
    Pass in either a `Species` or `Molecule` object and checks whether it passes
    the speciesCuttingThreshold set by the user. The default value is 20. If yes,
    returns `True` for passing cutting threshold.
    """

    from rmgpy.rmg.input import get_input

    try:
        species_constraints = get_input('species_constraints')
    except Exception:
        logging.debug('Species constraints could not be found.')
        species_constraints = {}

    if isinstance(species, Species):
        struct = species.molecule[0]
    else:
        # expects a molecule here
        struct = species

    min_cutting_size = species_constraints.get('speciesCuttingThreshold', 20)
    if struct.get_element_count()['C'] >= min_cutting_size:
        return True

    return False

def fails_species_constraints(species):
    """
    Pass in either a `Species` or `Molecule` object and checks whether it passes the
    speciesConstraints set by the user. If the species fails constraints, returns a string
    `reason` describing which constraint failed. If all constraints pass, returns `False`.

    Routing: polymer-sized species (or real Polymer objects) are bounded by the polymer tier;
    everything else by the global gas tier. See is_polymer_constraint_member.
    """
    from rmgpy.rmg.input import get_input

    try:
        polymer_constraints = get_input('polymer_constraints')
    except Exception:
        logging.debug('Polymer constraints could not be found.')
        polymer_constraints = None

    if is_polymer_constraint_member(species, polymer_constraints):
        if polymer_constraints is None:
            _warn_unbounded_polymer_once()
            return False
        struct = _normalize(species)
        reason = _evaluate_constraints(struct, polymer_constraints)
        if reason:
            _record_refusal(struct, 'polymer', reason)
        return reason

    try:
        species_constraints = get_input('species_constraints')
    except Exception:
        logging.debug('Species constraints could not be found.')
        species_constraints = {}

    struct = _normalize(species)
    reason = _evaluate_constraints(struct, species_constraints)
    if reason:
        _record_refusal(struct, 'gas', reason)
    return reason


def validate_explicit_dp_oligomers(initial_species, species_constraints):
    """
    Hard-error gate (never silent) for the explicit-DP handshake (stage A).

    Every auto-generated DP=xs oligomer (created by the polymer() input block
    when ``explicit_dp=True`` and marked with ``props['explicit_dp_origin']``,
    the sanctioned marker dict -- Species is a compiled extension type with no
    ad-hoc attributes) must
    survive the input-species constraint pass -- if the active species
    constraints would exclude it from the core, the handshake target would be
    silently absent and the feature structurally inert. Called from
    RMG.initialize BEFORE the generic input-species constraint loop so this
    actionable message (naming the deck flag) wins over the generic
    "remove the species" one, which is misleading for a species the user
    never wrote.

    ``'input species'`` in the constraints' ``allowed`` list admits the
    oligomer through the same escape hatch every other input species gets,
    so no error is raised in that case.
    """
    from rmgpy.exceptions import ForbiddenStructureException

    constraints = species_constraints or {}
    allowed = constraints.get('allowed', [])
    for spc in initial_species:
        origin = (getattr(spc, 'props', None) or {}).get('explicit_dp_origin')
        if origin is None:
            continue
        reason = _evaluate_constraints(_normalize(spc), constraints)
        if reason and 'input species' not in allowed:
            pool_label, dp = origin
            formula = spc.molecule[0].get_formula() if getattr(spc, 'molecule', None) else '?'
            raise ForbiddenStructureException(
                f"explicit_dp=True on polymer pool '{pool_label}' auto-generated "
                f"the DP={dp} capped oligomer '{spc.label}' ({formula}), but the "
                f"active species constraints refuse it: {reason}. The explicit-DP "
                f"handshake cannot run without this core species. Raise the "
                f"refusing maximum* limit in generatedSpeciesConstraints, add "
                f"'input species' to its 'allowed' list, or disable explicit_dp "
                f"on pool '{pool_label}'.")
