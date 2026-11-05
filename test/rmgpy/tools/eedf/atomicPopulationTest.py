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

"""Resolved atomic nitrogen composition axes and exact solver realization."""

import copy
import itertools

import numpy as np
import pytest

from rmgpy.tools.eedf.atomic_population import (
    AXES,
    BOUNDS,
    CONVENTION,
    STATES,
    validate_atomic_axes,
    realize_atomic_populations,
    resolve_atomic_species,
    atomic_held_out_points,
)
from rmgpy.tools.eedf.schema import SpecError, _content_hash


def document(refined=()):
    return {
        "axes": {
            name: ([0.0, top / 2, top] if name in refined else [0.0, top])
            for name, top in zip(AXES, BOUNDS)
        },
        "atomic_nitrogen_population": copy.deepcopy(CONVENTION),
    }


@pytest.mark.parametrize("refined", [(), AXES[:1], AXES])
def test_axis_policy_permits_only_one_midpoint_per_axis(refined):
    validate_atomic_axes(document(refined))


@pytest.mark.parametrize(
    "mutation",
    [
        "missing_policy",
        "missing_axis",
        "wrong_top",
        "extra_node",
        "nonmidpoint",
        "coordinate",
        "floor",
    ],
)
def test_atomic_axis_policy_refuses_incomplete_or_changed_convention(mutation):
    doc = document()
    if mutation == "missing_policy":
        del doc["atomic_nitrogen_population"]
    elif mutation == "missing_axis":
        del doc["axes"][AXES[1]]
    elif mutation == "wrong_top":
        doc["axes"][AXES[0]][-1] *= 2
    elif mutation == "extra_node":
        doc["axes"][AXES[0]] = [0.0, BOUNDS[0] / 4, BOUNDS[0] / 2, BOUNDS[0]]
    elif mutation == "nonmidpoint":
        doc["axes"][AXES[1]].insert(1, BOUNDS[1] / 3)
    elif mutation == "coordinate":
        doc["atomic_nitrogen_population"]["interpolation"] = "pchip"
    else:
        doc["atomic_nitrogen_population"]["zero_target_population_floor"] = 1e-5
    with pytest.raises(SpecError, match="atomic"):
        validate_atomic_axes(doc)


@pytest.mark.parametrize(
    "values", list(itertools.product(*[(0.0, top) for top in BOUNDS]))
)
def test_atomic_populations_preserve_total_gas_and_bounded_floor(values):
    gas = {"Ar": 0.95, "N2": 0.05}
    populations = {"Ar(1S0)": 1.0, "N2(X)": 1.0}
    actual_gas, actual_populations = realize_atomic_populations(
        document(), dict(zip(AXES, values)), gas, populations
    )
    assert gas == {"Ar": 0.95, "N2": 0.05}  # No caller mutation.
    assert actual_gas["N"] == sum(values)
    assert actual_gas["N2"] == 0.05
    assert sum(actual_gas.values()) == pytest.approx(1.0, abs=2e-16)
    assert sum(actual_populations[s] for s in STATES) == pytest.approx(
        1.0 if sum(values) else 0.0, abs=2e-16,
    )
    for state, requested in zip(STATES, values):
        assert actual_gas["N"] * actual_populations[state] == pytest.approx(
            requested, rel=1e-15, abs=1e-23
        )
    if any(values):
        # The pinned solver must see every forward target to assemble inverses.
        assert all(actual_populations[s] > 0 for s in STATES)


@pytest.mark.parametrize("value", [-1e-6, np.nan, np.inf, BOUNDS[0] * 1.00001])
def test_atomic_populations_refuse_outside_named_bound(value):
    with pytest.raises(SpecError, match=r"x_N\(4S\)"):
        realize_atomic_populations(
            document(), dict(zip(AXES, [value, 0.0, 0.0])), {"Ar": 0.95, "N2": 0.05}, {}
        )


def test_atomic_population_policy_is_part_of_row_identity():
    from rmgpy.tools.eedf.schema import _row_inputs

    base = {
        name: {}
        for name in (
            "arm",
            "Tg_K",
            "P_Pa",
            "loki_commit",
            "axes",
            "gas_properties",
            "state_properties",
            "working_conditions",
            "solver_options",
            "screen",
            "envelopes",
            "floors",
            "tolerances",
            "held_out",
            "refinement",
        )
    }
    base.update(document())
    inputs = _row_inputs(base)
    assert inputs["atomic_nitrogen_population"] == CONVENTION
    changed = copy.deepcopy(inputs)
    changed["atomic_nitrogen_population"]["zero_target_population_floor"] *= 10
    assert _content_hash(changed) != _content_hash(inputs)


def test_atomic_solver_realization_and_metadata_use_identical_populations():
    from rmgpy.tools.eedf.integrity import _resolved_properties

    spec = dict(
        document(),
        input_files={},
        solver_options={},
        envelopes={},
        gas_properties={"fraction": ["Ar = 0.95", "N2 = 0.05"]},
        state_properties={"population": ["Ar(1S0) = 1", "N2(X) = 1"]},
    )
    point = dict(zip(AXES, [0.0, BOUNDS[1], BOUNDS[2] / 2]))
    gas, props, fractions, pops = _resolved_properties(spec, point)
    assert {
        s: float(v) for s, v in (line.split(" = ") for line in gas["fraction"])
    } == fractions
    assert {
        s: float(v) for s, v in (line.split(" = ") for line in props["population"])
    } == pops
    assert fractions["N"] == BOUNDS[1] + BOUNDS[2] / 2
    assert fractions["N"] * pops["N(2D)"] == pytest.approx(BOUNDS[1], rel=1e-15)


def test_atomic_species_binding_uses_electronic_identity_not_labels():
    from rmgpy.species import Species

    s = Species(label="renamed S").from_adjacency_list("multiplicity 4\n1 N u3 p1 c0")
    d = Species(label="renamed D").from_adjacency_list(
        "multiplicity 2\nelectronicstate 2D\n1 N u1 p2 c0"
    )
    p = Species(label="renamed P").from_adjacency_list(
        "multiplicity 2\nelectronicstate 2P\n1 N u1 p2 c0"
    )
    assert resolve_atomic_species([p, s, d]) == dict(zip(AXES, [1, 2, 0]))
    with pytest.raises(SpecError, match="unique"):
        resolve_atomic_species([s, d])
    with pytest.raises(SpecError, match="unique"):
        resolve_atomic_species([s, d, p, copy.deepcopy(p)])


def test_atomic_held_out_covers_individual_and_joint_midpoints_and_refinement():
    axes = {
        "u": np.log([4.0, 17.0, 25.0]).tolist(),
        "x_N2": [0.01, 0.05, 0.2],
        "f_Ar4s": [0.0, 1e-5],
        **document()["axes"],
    }
    points = atomic_held_out_points(axes)
    # All 3^3-2^3 nonvertex midpoint configurations at every prescribed bath/field.
    assert len(points) == (27 - 8) * 3 * 3 * 2
    assert len({_content_hash(p) for p in points}) == len(points)
    for name, top in zip(AXES, BOUNDS):
        assert any(
            p[name] == top / 2 and all(p[n] == 0 for n in AXES if n != name)
            for p in points
        )
    assert any(all(p[n] == top / 2 for n, top in zip(AXES, BOUNDS)) for p in points)
    axes[AXES[0]] = [0.0, BOUNDS[0] / 2, BOUNDS[0]]
    refined = atomic_held_out_points(axes)
    assert {p[AXES[0]] for p in refined} >= {BOUNDS[0] / 4, 3 * BOUNDS[0] / 4}


def test_atomic_axes_refuse_generic_refinement_that_would_insert_thirds():
    doc = document()
    doc["refinement"] = {"max_rounds": 1}
    with pytest.raises(SpecError, match="T10-directed midpoint"):
        validate_atomic_axes(doc)


def test_screen_cannot_drop_declared_atomic_axes(monkeypatch):
    import rmgpy.tools.eedf.generate as generation

    doc = document()
    doc["axes"] = {"u": [1.0, 2.0], **doc["axes"]}
    doc["screen"] = {
        "threshold_fraction": 0.1,
        "candidates": [
            {"name": name, "reference": 0.0, "values": [0.0, top]}
            for name, top in zip(AXES, BOUNDS)
        ],
    }
    doc["envelopes"], doc["tolerances"], doc["floors"] = {}, {}, {}
    doc["Tg_K"] = 300.0
    monkeypatch.setattr(
        generation, "_solve", lambda *args: [{"u": u} for u in [1.0, 2.0]]
    )
    monkeypatch.setattr(generation, "_compare_row", lambda *args: [])
    results = generation._screen(doc, None, [], "fixture")
    assert all(result["selected_axis"] for result in results)
    assert set(AXES).issubset(doc["axes"])


def test_standard_held_out_selection_includes_atomic_midpoints():
    from rmgpy.tools.eedf.validation import _held_out_points

    axes = {"u": [1.0, 2.0], **document()["axes"]}
    actual = _held_out_points(
        axes, {"cell_midpoints": False, "lhs_count": 0, "seed": 17}
    )
    assert {_content_hash(p) for p in atomic_held_out_points(axes)} <= {
        _content_hash(p) for p in actual
    }


def test_atomic_t10_requires_frozen_accuracy_and_complete_midpoint_coverage():
    from rmgpy.tools.eedf.atomic_population import validate_atomic_qualification

    doc = document()
    doc["axes"] = {"u": [1.0, 2.0], **doc["axes"]}
    criteria = {
        name: {"rtol": 0.01, "atol": 0.0} for name in ("mobility_N", "field_power")
    }
    points = atomic_held_out_points(doc["axes"])
    validate_atomic_qualification(doc, criteria, points)
    with pytest.raises(SpecError, match="atomic.*criteria"):
        validate_atomic_qualification(
            doc, {name: {"rtol": 0.02, "atol": 0.0} for name in criteria}, points
        )
    with pytest.raises(SpecError, match="atomic.*midpoint"):
        validate_atomic_qualification(doc, criteria, points[:-1])


def test_atomic_refinement_changes_only_requested_axis_and_requires_fresh_midpoints():
    from rmgpy.tools.eedf.atomic_population import refine_atomic_axes

    doc = document()
    doc["axes"] = {"u": [1.0, 2.0], **doc["axes"]}
    doc["held_out"] = {
        "cell_midpoints": True,
        "lhs_count": 0,
        "seed": 1,
        "off_centre": atomic_held_out_points(doc["axes"]),
    }
    revised = refine_atomic_axes(doc, [AXES[0]])
    assert doc["axes"][AXES[0]] == [0.0, BOUNDS[0]]
    assert revised["axes"][AXES[0]] == [0.0, BOUNDS[0] / 2, BOUNDS[0]]
    assert revised["axes"][AXES[1]] == doc["axes"][AXES[1]]
    assert revised["axes"][AXES[2]] == doc["axes"][AXES[2]]
    native = set(itertools.product(*revised["axes"].values()))
    assert all(
        tuple(p[name] for name in revised["axes"]) not in native
        for p in revised["held_out"]["off_centre"]
    )
    with pytest.raises(SpecError, match="already.*refined"):
        refine_atomic_axes(revised, [AXES[0]])
