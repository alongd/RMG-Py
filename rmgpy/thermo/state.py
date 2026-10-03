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

"""Library-only thermochemistry for resolved states and declared v=0 species."""

from rmgpy.exceptions import ExcitedSpeciesThermoError, VibrationalManifoldError


def require_thermo_estimation_allowed(species_or_molecule):
    """Refuse state-blind estimation before a caller mutates or converts a graph."""
    from rmgpy.species import TransitionState
    if isinstance(species_or_molecule, TransitionState):
        return
    molecules = getattr(species_or_molecule, 'molecule', [species_or_molecule])
    label = getattr(species_or_molecule, 'label', '')
    if species_or_molecule.props.get('vibrational_manifold'):
        raise VibrationalManifoldError(
            "Species {0!r} is declared by vibrationalManifold and requires a loaded "
            "thermo library entry for the same graph with vibrationallevel 0; "
            "other thermo sources are refused.".format(label or species_or_molecule.props['vibrational_manifold']))
    for molecule in molecules:
        if molecule.props.get('vibrational_manifold'):
            raise VibrationalManifoldError(
                "Species {0!r} is declared by vibrationalManifold; thermo must come "
                "from a library entry with vibrationallevel 0.".format(label or molecule.props['vibrational_manifold']))
        if molecule.has_resolved_state():
            raise ExcitedSpeciesThermoError(
                "Library-only thermochemistry is required for resolved-state species {0!r}. "
                "A loaded thermo library must contain an entry for this exact state; "
                "QM, HBI, ML and group additivity cannot estimate its thermo.\n{1}"
                .format(label or molecule.to_adjacency_list().strip(), molecule.to_adjacency_list()))


def require_consistent_state(species):
    """Refuse inconsistent public state headers across resonance candidates."""
    states = {(mol.electronic_state, mol.vibrational_level) for mol in species.molecule}
    if len(states) > 1:
        raise ExcitedSpeciesThermoError(
            "Species {0!r} has inconsistent resolved states across its structures; "
            "thermo requires one exact library state.".format(species.label))
    if any(mol.props.get('vibrational_manifold') for mol in species.molecule) and not species.props.get('vibrational_manifold'):
        raise VibrationalManifoldError(
            "Species {0!r} has a molecule-level vibrationalManifold declaration without "
            "a Species declaration; use vibrationalManifold before thermo lookup.".format(species.label))


def thermo_library_species(species):
    """Return the library lookup structure without changing model identity."""
    require_consistent_state(species)
    if not species.props.get('vibrational_manifold'):
        return species
    lookup = species.copy(deep=True)
    lookup.props.pop('vibrational_manifold', None)
    for molecule in lookup.molecule:
        molecule.props.pop('vibrational_manifold', None)
        molecule.vibrational_level = 0
    return lookup


def require_electron_state_allowed(species_or_molecule):
    """The charge-carrier electron cannot carry molecular state keywords."""
    molecules = getattr(species_or_molecule, 'molecule', [species_or_molecule])
    for molecule in molecules:
        if molecule.is_electron() and (molecule.has_resolved_state()
                or species_or_molecule.props.get('vibrational_manifold')
                or molecule.props.get('vibrational_manifold')):
            raise ExcitedSpeciesThermoError(
                "Species {0!r}: electron electronicstate, vibrationallevel and vibrationalManifold "
                "are refused; the electron requires its canonical unkeyed library entry.\n{1}"
                .format(getattr(species_or_molecule, 'label', 'electron'), molecule.to_adjacency_list()))


def requires_state_library(species):
    """Whether attached or generated thermo needs exact-state library provenance."""
    require_consistent_state(species)
    return (bool(species.props.get('vibrational_manifold'))
            or any(molecule.has_resolved_state() for molecule in species.molecule))


def _thermo_numbers_match(actual, expected):
    """Compare finite structural coefficients, including exact zero coefficients."""
    from math import isclose, isfinite
    return (isfinite(actual) and isfinite(expected)
            and isclose(actual, expected, rel_tol=1e-9, abs_tol=0.0))


def _thermo_quantity_match(actual, expected):
    """Compare optional SI quantities without treating missing fields as zero."""
    if actual is None or expected is None:
        return actual is expected
    return _thermo_numbers_match(actual.value_si, expected.value_si)


def thermo_fields_match(thermo, reference):
    """Compare inspectable thermo structures; samples never authorize a source."""
    from rmgpy.thermo import NASA, NASAPolynomial, ThermoData, Wilhoit
    if thermo is None or reference is None:
        return False
    supported = (NASA, NASAPolynomial, ThermoData, Wilhoit)
    if type(thermo) not in supported or type(reference) not in supported:
        raise ExcitedSpeciesThermoError(
            "Cannot inspect thermo model classes {!r} and {!r} for exact library-source equality."
            .format(type(thermo).__name__, type(reference).__name__))
    if type(thermo) is not type(reference):
        return False
    if getattr(thermo, 'thermo_coverage_dependence', {}) != getattr(reference, 'thermo_coverage_dependence', {}):
        return False
    for name in ('E0', 'Cp0', 'CpInf', 'Tmin', 'Tmax'):
        if not _thermo_quantity_match(getattr(thermo, name, None), getattr(reference, name, None)):
            return False
    if type(thermo) is NASA:
        from rmgpy.solver.plasma import _thermo_coverage_gap
        if (_thermo_coverage_gap(thermo) is not None
                or _thermo_coverage_gap(reference) is not None
                or not thermo.polynomials
                or len(thermo.polynomials) != len(reference.polynomials)):
            return False
        return all(thermo_fields_match(actual, expected)
                   for actual, expected in zip(thermo.polynomials, reference.polynomials))
    if type(thermo) is NASAPolynomial:
        actual, expected = thermo.coeffs, reference.coeffs
        return (len(actual) == len(expected)
                and all(_thermo_numbers_match(a, e) for a, e in zip(actual, expected)))
    if type(thermo) is ThermoData:
        for name in ('Tdata', 'Cpdata'):
            actual, expected = getattr(thermo, name).value_si, getattr(reference, name).value_si
            if len(actual) != len(expected) or not all(
                    _thermo_numbers_match(a, e) for a, e in zip(actual, expected)):
                return False
        return all(_thermo_quantity_match(getattr(thermo, name), getattr(reference, name))
                   for name in ('H298', 'S298'))
    return (all(_thermo_numbers_match(getattr(thermo, name), getattr(reference, name))
                for name in ('a0', 'a1', 'a2', 'a3'))
            and all(_thermo_quantity_match(getattr(thermo, name), getattr(reference, name))
                    for name in ('B', 'H0', 'S0')))


def require_state_thermo(species, thermo_database, solvent_name='', thermo=None):
    """Value-check attached thermo against exact-state library data and processing.

    Cached values carry no authority: a matching entry must still be loaded.
    Derived solvation and state-blind limits are refused.
    """
    if not requires_state_library(species):
        return
    require_electron_state_allowed(species)
    if solvent_name:
        require_thermo_estimation_allowed(species)
    if thermo is None:
        thermo = species.thermo
    from copy import deepcopy
    from rmgpy.data.thermo import find_cp0_and_cpinf
    from rmgpy.thermo.thermoengine import _process_thermo_data
    error = (VibrationalManifoldError if species.props.get('vibrational_manifold')
             else ExcitedSpeciesThermoError)
    lookup = thermo_library_species(species)
    if thermo_database is not None:
        for name in thermo_database.library_order:
            library = thermo_database.libraries[name]
            if library.solvent and library.solvent != solvent_name:
                continue
            for entry in library.entries.values():
                if entry.data is None or not any(mol.is_isomorphic(entry.item) for mol in lookup.molecule):
                    continue
                reference = deepcopy(entry.data)
                reference.label = entry.label
                find_cp0_and_cpinf(lookup, reference)
                reference.comment += ('Liquid thermo library: ' if library.solvent else 'Thermo library: ') + name
                forms = [reference]
                forms.append(_process_thermo_data(lookup.copy(deep=True), deepcopy(reference),
                                                 solvent_name=solvent_name))
                for form in forms:
                    if thermo_fields_match(thermo, form):
                        return form
    raise error(
        "Thermo for species {0!r} does not value-match a loaded thermo library entry "
        "for its exact state{1}; attached or cached thermo cannot bypass library-only "
        "thermochemistry.\n{2}".format(
            species.label, ' with vibrationallevel 0' if species.props.get('vibrational_manifold') else '',
            lookup.molecule[0].to_adjacency_list()))


def require_vibrational_manifold_thermo(species, thermo_database):
    """Check attached v=0 thermo against its matching loaded library entry."""
    require_state_thermo(species, thermo_database)


def _needs_state_thermo_check(species):
    """Include formerly resolved caches until the getter invalidates their source."""
    from rmgpy.species import TransitionState
    if isinstance(species, TransitionState):
        return False
    if requires_state_library(species):
        return True
    previous = species.props.get('_thermo_state_key')
    return (previous is not None
            and (bool(previous[1]) or any(es or v >= 0 for es, v in previous[0])))


def checked_thermo(species):
    """Check current or former resolved owners; preserve ordinary ground attachments."""
    return species.get_thermo_data() if _needs_state_thermo_check(species) else species.thermo


def checked_energy(species):
    """Refresh current or former resolved energy from the checked thermo source."""
    if _needs_state_thermo_check(species):
        species.set_e0_with_thermo()
    return species.conformer.E0.value_si


def require_species_thermo_allowed(species_list):
    """Refuse unsupported derivation or statmech for every participating species."""
    for species in species_list:
        require_thermo_estimation_allowed(species)


def require_network_thermo_allowed(network):
    """Pressure-dependence statmech is unsupported for library-only states."""
    for configuration in network.isomers + network.reactants + network.products:
        require_species_thermo_allowed(configuration.species)
    for reaction in network.path_reactions:
        require_species_thermo_allowed(reaction.reactants + reaction.products)


# Atom-list heuristics have no state header. Retain weak ownership so extracting
# a cycle cannot silently erase its molecule's resolved identity. No state is
# stored on atoms, and shared atoms are checked against every live owner.
_resolved_atom_owners = {}


def register_state_atoms(molecule):
    """Register live resolved graph ownership without changing atom data or persistence."""
    if not (molecule.has_resolved_state() or molecule.props.get('vibrational_manifold')):
        return
    import weakref
    owner_id = id(molecule)
    atom_ids = tuple(id(atom) for atom in molecule.vertices)
    def discard(reference):
        for atom_id in atom_ids:
            owners = _resolved_atom_owners.get(atom_id)
            if owners is not None and owners.get(owner_id) is reference:
                del owners[owner_id]
                if not owners:
                    del _resolved_atom_owners[atom_id]
    reference = weakref.ref(molecule, discard)
    for atom_id in atom_ids:
        _resolved_atom_owners.setdefault(atom_id, {})[owner_id] = reference


def require_atom_thermo_allowed(atoms):
    """Refuse atom-only estimation if any atom belongs to a resolved graph."""
    for atom in atoms:
        for reference in tuple(_resolved_atom_owners.get(id(atom), {}).values()):
            owner = reference()
            if owner is not None and any(vertex is atom for vertex in owner.vertices):
                require_thermo_estimation_allowed(owner)
