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

"""Resolved neutral nitrogen fractions for the atomic collision targets."""

import copy
import itertools

import numpy as np

AXES = ("x_N(4S)", "x_N(2D)", "x_N(2P)")
STATES = ("N(4S)", "N(2D)", "N(2P)")
BOUNDS = (5.0e-4, 2.0e-4, 2.0e-4)
CONVENTION = {
    "axes": list(AXES),
    "states": list(STATES),
    "fraction_basis": "total-neutral-gas",
    "displaced_gas": "Ar",
    "interpolation": "multilinear",
    "zero_target_population_floor": 1.0e-20,
    "zero_gas_populations": "all-zero",
}


def validate_atomic_axes(document):
    """Require all three axes and the frozen physical/interpolation convention.

    Each axis may contain its two endpoints or those endpoints plus one
    midpoint. A further refinement requires a new scientific decision.
    """
    from rmgpy.tools.eedf.schema import SpecError

    axes = document.get("axes", {})
    declaration = document.get("atomic_nitrogen_population")
    if declaration is None and not set(AXES).intersection(axes):
        return
    if declaration != CONVENTION or not set(AXES).issubset(axes):
        raise SpecError(
            "atomic nitrogen axes require their complete declared convention"
        )
    if document.get("refinement", {}).get("max_rounds", 0) != 0:
        raise SpecError("atomic axes require explicit T10-directed midpoint refinement")
    for name, bound in zip(AXES, BOUNDS):
        nodes = np.asarray(axes[name], dtype=float)
        if not (
            np.array_equal(nodes, [0.0, bound])
            or np.array_equal(nodes, [0.0, bound / 2.0, bound])
        ):
            raise SpecError(
                "atomic axis "
                + name
                + " permits only endpoints and one midpoint refinement"
            )


def realize_atomic_populations(document, coordinates, gas, populations):
    """Realize total-gas atomic fractions, displacing the same amount of Ar.

    The within-N population floor avoids the pinned LoKI inverse-assembly
    defect at a zero forward target with a populated inverse source. Its
    maximum total-gas perturbation is below 1e-23. The identical floored
    populations are used for solver input and energy-schema row metadata.
    With no atomic gas all electronic populations and collision densities
    remain zero: LoKI gives the absent gas a zero root population and refuses
    populated children under that root.
    """
    from rmgpy.tools.eedf.schema import SpecError

    validate_atomic_axes(document)
    values = []
    for name, bound in zip(AXES, BOUNDS):
        try:
            value = float(coordinates[name])
        except (KeyError, TypeError, ValueError) as error:
            raise SpecError("atomic coordinate missing or invalid: " + name) from error
        if not np.isfinite(value) or not 0.0 <= value <= bound:
            raise SpecError(
                "atomic axis {0} outside zero to top node {1:g}".format(name, bound)
            )
        values.append(value)
    total = sum(values)
    fractions = dict(gas)
    if "Ar" not in fractions:
        raise SpecError("atomic nitrogen population requires its displaced Ar gas")
    fractions["Ar"] += fractions.get("N", 0.0) - total
    fractions["N"] = total
    if fractions["Ar"] < 0.0 or not np.isclose(
        sum(fractions.values()), 1.0, rtol=0.0, atol=1.0e-12
    ):
        raise SpecError("atomic gas fractions must be normalized with nonnegative Ar")
    result = dict(populations)
    if total:
        conditional = np.maximum(
            np.asarray(values) / total, CONVENTION["zero_target_population_floor"]
        )
        conditional /= conditional.sum()
    else:
        conditional = (0.0, 0.0, 0.0)
    result.update(zip(STATES, map(float, conditional)))
    return fractions, result


def resolve_atomic_species(species):
    """Bind one neutral atom for each term by molecular electronic identity."""
    from rmgpy.tools.eedf.schema import SpecError

    found = {name: [] for name in AXES}
    identities = {("", 4): AXES[0], ("2D", 2): AXES[1], ("2P", 2): AXES[2]}
    for index, spc in enumerate(species):
        if spc.get_net_charge() or len(spc.molecule) != 1:
            continue
        molecule = spc.molecule[0]
        if len(molecule.atoms) != 1 or molecule.get_formula() != "N":
            continue
        key = (molecule.electronic_state, molecule.multiplicity)
        if key in identities:
            found[identities[key]].append(index)
    if any(len(indices) != 1 for indices in found.values()):
        raise SpecError(
            "atomic axes require a unique neutral N(4S), N(2D), N(2P) core species"
        )
    return {name: indices[0] for name, indices in found.items()}


def atomic_held_out_points(axes):
    """Return individual and joint atomic cell-midpoints at every other node.

    Include partially joint midpoints too, so cancellation at the all-midpoint
    state cannot hide curvature on a pair of atomic axes.
    """
    names = tuple(axes)
    if not set(AXES).issubset(names):
        from rmgpy.tools.eedf.schema import SpecError

        raise SpecError("atomic held-out points require all three axes")
    options = []
    for name in names:
        nodes = np.asarray(axes[name], dtype=float)
        options.append(
            sorted(
                set(nodes)
                | (set((nodes[:-1] + nodes[1:]) / 2.0) if name in AXES else set())
            )
        )
    native = set(itertools.product(*(axes[name] for name in names)))
    return [
        dict(zip(names, map(float, values)))
        for values in itertools.product(*options)
        if values not in native
    ]


def validate_atomic_qualification(document, criteria, points=None):
    """Require the ruled T10 tolerances and every atomic midpoint reference."""
    from rmgpy.tools.eedf.schema import SpecError

    if "atomic_nitrogen_population" not in document:
        return
    validate_atomic_axes(document)
    expected = {
        name: {"rtol": 0.01, "atol": 0.0} for name in ("mobility_N", "field_power")
    }
    if criteria != expected:
        raise SpecError("atomic T10 criteria require rtol=0.01 and atol=0")
    if points is not None:
        names = tuple(document["axes"])
        required = {
            tuple(point[name] for name in names)
            for point in atomic_held_out_points(document["axes"])
        }
        retained = {tuple(point[name] for name in names) for point in points}
        if not required.issubset(retained):
            raise SpecError("atomic T10 midpoint reference coverage is incomplete")


def refine_atomic_axes(document, failed_axes):
    """Insert one midpoint on each explicitly identified failing atomic axis.

    Return a new spec and fresh held-out coordinates; callers regenerate and
    requalify rather than promote a previously held-out reference silently.
    """
    from rmgpy.tools.eedf.schema import SpecError

    validate_atomic_axes(document)
    if not failed_axes or not set(failed_axes).issubset(AXES):
        raise SpecError("atomic refinement requires the named failing axes")
    result = copy.deepcopy(document)
    for name in failed_axes:
        axis = result["axes"][name]
        if len(axis) != 2:
            raise SpecError(
                "atomic axis " + name + " already midpoint-refined; needs input"
            )
        result["axes"][name] = [axis[0], (axis[0] + axis[-1]) / 2.0, axis[-1]]
    names = tuple(result["axes"])
    native = set(itertools.product(*(result["axes"][name] for name in names)))
    retained = {
        tuple(point[name] for name in names)
        for point in result.get("held_out", {}).get("off_centre", [])
    }
    retained.update(
        tuple(point[name] for name in names)
        for point in atomic_held_out_points(result["axes"])
    )
    if "held_out" in result:
        result["held_out"]["off_centre"] = [
            dict(zip(names, point)) for point in sorted(retained - native)
        ]
    return result
