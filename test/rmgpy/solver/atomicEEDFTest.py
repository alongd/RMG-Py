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

"""Atomic EEDF interpolation and reactor fractions, using actual table artifacts."""

import copy
import json
from types import SimpleNamespace

import rmgpy.data.rmg as rmg_data_module
from rmgpy.data.base import Entry
from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary

import numpy as np
import pytest

from rmgpy import constants
from rmgpy.solver.eedf import EEDFTable
from rmgpy.solver.plasma import PlasmaReactor, PlasmaStateError
from rmgpy.species import Species
from rmgpy.tools.eedf.artifact import _write_artifact
from rmgpy.tools.eedf.atomic_population import AXES, BOUNDS, CONVENTION
from rmgpy.tools.eedf.schema import _content_hash, _file_hash, FingerprintMismatch
from plasmaEEDFCoreTest import (
    artifact,
    species_and_reaction,
    thermo,
    TGAS,
    PRESSURE,
    VOLUME,
    LAMBDA,
    MU0,
)


@pytest.fixture(autouse=True)
def isolated_database(monkeypatch):
    monkeypatch.setattr(rmg_data_module, "database", None)


def atomic_artifact(tmp_path, marker, argon=False):
    seed_path = artifact(
        tmp_path / "seed",
        marker,
        composition=False,
        gases=["Ar", "N"],
        include_ar4s_envelope=False,
    )
    manifest = json.loads((seed_path / "manifest.json").read_text())
    seed = EEDFTable.load(
        seed_path,
        manifest["row_inputs"],
        artifact_sha256=_file_hash(seed_path / "table.h5"),
    )
    axes = {
        "u": [0.0, 1.0, 2.0, 3.0],
        **{
            name: ([0.0, top / 2, top] if i == 0 else [0.0, top])
            for i, (name, top) in enumerate(zip(AXES, BOUNDS))
        },
    }
    if argon:
        from rmgpy.tools.eedf.argon_population import NODES
        axes['f_Ar4s'] = list(NODES)
    axes = {'u': axes['u'], **{name: axes[name] for name in sorted(axes) if name != 'u'}}
    rows = np.empty(tuple(len(v) for v in axes.values()), dtype=object)
    for index in np.ndindex(rows.shape):
        u = axes["u"][index[0]]
        point = {
            name: axes[name][index[i + 1]]
            for i, name in enumerate(name for name in axes if name != "u")
        }
        s, d, p = [point[name] / top for name, top in zip(AXES, BOUNDS)]
        row = seed.row(u, {})._as_dict()
        row["composition"] = point
        row["k_ine"] *= 1 + s * s
        row["swarm"]["mean_energy_eV"] = 1 + u + s * s + 2 * d + 3 * p + 0.4 * s * d
        if argon:
            from rmgpy.tools.eedf.argon_population import coordinate
            a = float(coordinate(point['f_Ar4s']))
            row['swarm']['mean_energy_eV'] += (0.2 + 0.1 * d) * a
        row["swarm"]["mobility_N"] = 1e24 * (
            1 + 0.15 * u + 0.1 * s * s + 0.2 * d + 0.3 * p + 0.04 * s * d
        )
        if argon:
            row['swarm']['mobility_N'] += 1e24 * (0.02 + 0.01*d) * a
        row["power_groups"]["field"] = (
            row["swarm"]["mobility_N"] * (np.exp(u) * 1e-21) ** 2
        )
        rows[index] = row
    model = copy.deepcopy(manifest["row_inputs"])
    model["axes"] = axes
    model["atomic_nitrogen_population"] = copy.deepcopy(CONVENTION)
    if argon:
        from rmgpy.tools.eedf.argon_population import CONVENTION as argon_convention
        model['argon_4s_population'] = copy.deepcopy(argon_convention)
    manifest.update(
        axes=axes,
        row_inputs=model,
        fingerprint=_content_hash(model),
        held_out=[],
        held_out_verdicts=[
            {
                "point": {"u": 0.5},
                "branch_id": "branch_0",
                "passed": True,
                "checks": [{"passed": True}],
            }
        ],
    )
    path = _write_artifact(tmp_path / "atomic", manifest, {"branch_0": rows}, [])
    return path, model


@pytest.mark.parametrize('fraction', [3e-8, 3e-6])
def test_atomic_and_argon_axes_share_interpolation_and_chain_rule(tmp_path, fraction):
    from rmgpy.tools.eedf.argon_population import coordinate

    _, reaction = species_and_reaction()
    path, model = atomic_artifact(tmp_path, reaction.kinetics, argon=True)
    table = EEDFTable.load(path, model, artifact_sha256=_file_hash(path / 'table.h5'))
    s, d, p, u = 0.25, 0.3, 0.6, 0.7
    point = {name: top * value for name, top, value in zip(AXES, BOUNDS, (s, d, p))}
    point['f_Ar4s'] = fraction
    row = table.row(u, point)
    a = float(coordinate(fraction))
    expected = 1 + u + 0.125 + 2*d + 3*p + 0.4*s*d + (0.2 + 0.1*d)*a
    assert row.swarm['mean_energy_eV'] == pytest.approx(expected, rel=1e-13)
    derivatives = row.derivatives['composition']
    assert derivatives['f_Ar4s']['swarm']['mobility_N'] == pytest.approx(
        1e24 * (0.02 + 0.01*d) / max(fraction, 1e-7), rel=1e-12,
    )
    assert derivatives[AXES[1]]['swarm']['mobility_N'] == pytest.approx(
        1e24 * (0.2 + 0.04*s + 0.01*a) / BOUNDS[1], rel=1e-12,
    )


def test_atomic_table_is_multilinear_even_after_midpoint_refinement(tmp_path):
    _, reaction = species_and_reaction()
    path, model = atomic_artifact(tmp_path, reaction.kinetics)
    table = EEDFTable.load(path, model, artifact_sha256=_file_hash(path / "table.h5"))
    s, d, p, u = 0.25, 0.3, 0.6, 0.7
    point = {name: top * v for name, top, v in zip(AXES, BOUNDS, (s, d, p))}
    row = table.row(u, point)
    # Piecewise linear s^2 between nodes 0 and .5 gives .125, not PCHIP.
    expected_energy = 1 + u + 0.125 + 2 * d + 3 * p + 0.4 * s * d
    expected_mobility = 1e24 * (
        1 + 0.15 * u + 0.1 * 0.125 + 0.2 * d + 0.3 * p + 0.04 * s * d
    )
    assert row.swarm["mean_energy_eV"] == pytest.approx(expected_energy, rel=1e-14)
    assert row.swarm["mobility_N"] == pytest.approx(expected_mobility, rel=1e-14)
    assert row.power_groups["field"] == pytest.approx(
        expected_mobility * (np.exp(u) * 1e-21) ** 2
    )
    rate_row = table.row(1.0, point)
    left, right = dict(point), dict(point)
    left[AXES[0]], right[AXES[0]] = 0.0, BOUNDS[0] / 2
    expected_rates = (table.row(1.0, left).k_ine + table.row(1.0, right).k_ine) / 2
    assert rate_row.k_ine == pytest.approx(expected_rates, rel=1e-14)
    for name, bound in zip(AXES, BOUNDS):
        delta = bound * 1e-5
        low, high = dict(point), dict(point)
        low[name] -= delta
        high[name] += delta
        fd = (
            table.row(u, high).power_groups["field"]
            - table.row(u, low).power_groups["field"]
        ) / (2 * delta)
        assert row.derivatives["composition"][name]["power_groups"][
            "field"
        ] == pytest.approx(fd, rel=2e-9)


def test_atomic_table_refuses_axis_without_population_policy(tmp_path):
    _, reaction = species_and_reaction()
    path, model = atomic_artifact(tmp_path, reaction.kinetics)
    manifest = json.loads((path / "manifest.json").read_text())
    del manifest["row_inputs"]["atomic_nitrogen_population"]
    manifest["fingerprint"] = _content_hash(manifest["row_inputs"])
    (path / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises((FingerprintMismatch, ValueError), match="atomic"):
        EEDFTable.load(
            path, manifest["row_inputs"], artifact_sha256=_file_hash(path / "table.h5")
        )


def atomic_reactor(tmp_path):
    core, reaction = species_and_reaction()
    for label, adjacency in [
        ("renamed ground", "multiplicity 4\n1 N u3 p1 c0"),
        ("renamed doublet D", "multiplicity 2\nelectronicstate 2D\n1 N u1 p2 c0"),
        ("renamed doublet P", "multiplicity 2\nelectronicstate 2P\n1 N u1 p2 c0"),
    ]:
        spc = Species(label=label).from_adjacency_list(adjacency)
        spc.thermo = thermo()
        core.append(spc)
    database = ThermoDatabase()
    atomic_library = ThermoLibrary(label="SyntheticAtomic")
    argon_library = ThermoLibrary(label="SyntheticArgon", thermo_convention="ion")
    for index, species in enumerate(core[1:], 1):
        species.thermo.Cp0 = (20.0, "J/(mol*K)")
        species.thermo.CpInf = (20.0, "J/(mol*K)")
        library = atomic_library if index >= 4 else argon_library
        library.entries[index] = Entry(
            index=index,
            label=species.label,
            item=species.molecule[0].copy(deep=True),
            data=copy.deepcopy(species.thermo),
        )
    database.libraries = {lib.label: lib for lib in (atomic_library, argon_library)}
    database.library_order = [atomic_library.label, argon_library.label]
    rmg_data_module.database = SimpleNamespace(thermo=database, solvation=None)
    path, model = atomic_artifact(tmp_path, reaction.kinetics)
    config = {
        "provider": "loki-table",
        "table": (str(path), _file_hash(path / "table.h5")),
        "branch": "branch_0",
        "initial_reduced_field": (np.e, "Td"),
        "empirical_laws": {},
    }
    initial = {spc: 0.0 for spc in core}
    initial.update({core[0]: 1e-6, core[1]: 0.999998, core[2]: 1e-6})
    reactor = PlasmaReactor(
        (TGAS, "K"),
        (PRESSURE, "Pa"),
        initial,
        (TGAS, "K"),
        termination=[],
        diffusion_length=(LAMBDA, "m"),
        ion_reduced_mobility=(MU0, "m^2/(V*s)"),
        wall_single_bath_approximation=True,
        wall_neutralization_products={"Ar+": "Ar"},
        electron_energy_balance={
            "absorbed_power": (0.5, "W"),
            "chamber_volume": (VOLUME, "m^3"),
            "sheath": "floating_wall",
        },
        electron_kinetics=config,
        thermo_source_assertions={"Ar+": "ion"},
    )
    reactor.initialize_model(core, [reaction], [], [])
    return reactor, core


def atomic_trial(reactor, values):
    y = reactor.y.copy()
    neutral = y[1] + sum(y[3 : reactor.num_core_species])
    y[1] -= sum(values) * neutral
    y[4:7] = np.asarray(values) * neutral
    return y


def test_atomic_reactor_realizes_renamed_species_and_energy_chain_rule(tmp_path):
    reactor, core = atomic_reactor(tmp_path)
    values = np.asarray(BOUNDS) * [0.25, 0.3, 0.6]
    y = atomic_trial(reactor, values)
    point = reactor._eedf_coordinates(y)
    assert [point[name] for name in AXES] == pytest.approx(values, rel=1e-14)
    rates = np.zeros_like(y)
    rates[4:7] = [2e-8, -1e-8, 3e-8]
    calculated = reactor._eedf_composition_energy_derivative(y, rates)
    step = 1e-3
    provider = reactor.eedf_provider
    energies = [
        provider.mean_energy_eV(
            provider.row(
                float(y[reactor.te_index]),
                reactor._eedf_coordinates(y + sign * step * rates),
            )
        )
        for sign in (-1, 1)
    ]
    finite_difference = (
        y[reactor.electron_index]
        * constants.Na
        * constants.e
        * (energies[1] - energies[0])
        / (2 * step)
    )
    assert calculated == pytest.approx(finite_difference, rel=2e-6, abs=1e-13)


@pytest.mark.parametrize("index", range(3))
@pytest.mark.parametrize("clamp", [False, True])
def test_atomic_reactor_refuses_named_top_node_even_on_trial_path(
    tmp_path, index, clamp
):
    reactor, _ = atomic_reactor(tmp_path)
    values = np.zeros(3)
    values[index] = BOUNDS[index] * 1.001
    with pytest.raises(PlasmaStateError) as exc:
        reactor._eedf_coordinates(atomic_trial(reactor, values), clamp=clamp)
    assert AXES[index] in str(exc.value)
    assert "top node" in str(exc.value)


@pytest.mark.parametrize("near_top", [False, True])
def test_atomic_coordinate_energy_jacobian_matches_residual(tmp_path, near_top):
    reactor, _ = atomic_reactor(tmp_path)
    values = np.asarray(BOUNDS) * ([1.0 - 1e-8] * 3 if near_top else [0.25, 0.3, 0.6])
    y = atomic_trial(reactor, values)
    dydt = np.zeros_like(y)
    baseline = reactor.residual(0.0, y, dydt)[0]
    jacobian = np.asarray(reactor.jacobian(0.0, y, dydt, 0.0))
    # The Ar denominator column also moves every atomic fraction. At the
    # upper face its admissible direction is opposite to an atomic column's.
    for column in (1, 4, 5, 6):
        step = 5e-6 * abs(y[column])
        below = y.copy()
        below[column] -= step
        above = y.copy()
        above[column] += step
        if near_top:
            if column == 1:
                finite_difference = (
                    reactor.residual(0.0, above, dydt)[0] - baseline
                ) / step
            else:
                finite_difference = (
                    baseline - reactor.residual(0.0, below, dydt)[0]
                ) / step
        else:
            finite_difference = (
                reactor.residual(0.0, above, dydt)[0]
                - reactor.residual(0.0, below, dydt)[0]
            ) / (2 * step)
        assert finite_difference[reactor.te_index] != 0.0
        scale = np.max(np.abs(finite_difference))
        np.testing.assert_allclose(
            jacobian[:, column], finite_difference, rtol=3e-4, atol=2e-6 * scale
        )


@pytest.mark.parametrize("index", range(3))
def test_atomic_jacobian_refuses_out_of_domain_base_state(tmp_path, index):
    reactor, _ = atomic_reactor(tmp_path)
    values = np.asarray(BOUNDS) * 0.5
    values[index] = BOUNDS[index] * 1.001
    y = atomic_trial(reactor, values)
    with pytest.raises(PlasmaStateError, match="top node") as exc:
        reactor.jacobian(0.0, y, np.zeros_like(y), 0.0)
    assert AXES[index] in str(exc.value)
