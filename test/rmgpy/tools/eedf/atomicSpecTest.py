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

"""Atomic source selection, ordinary owner binding, and spec realization."""

import copy
import hashlib
import json
import os
from pathlib import Path

import pytest

from rmgpy.tools.eedf.atomic import realize_sources, refresh_spec
from rmgpy.tools.eedf.channels import _classic_sections, _load_channel_map
from rmgpy.tools.eedf.schema import FingerprintMismatch, SpecError, load_spec
from rmgpy.tools.eedf.atomic_population import AXES, BOUNDS, CONVENTION


@pytest.fixture
def sources():
    source = os.environ.get("EEDF_ATOMIC_SOURCE_DIR")
    database = os.environ.get("EEDF_DATABASE_ROOT")
    if not source or not database:
        pytest.skip("pinned atomic source and database paths are required")
    return Path(source), Path(database)


def test_realization_retains_exact_source_rows_and_excludes_ionization(tmp_path, sources):
    source, database = sources
    manifest = database / "input/eedf/atomic_nitrogen_channel_map.json"
    result = realize_sources(source, manifest, tmp_path / "realization")
    selected = _classic_sections(result["path"], result["sha256"])
    original = []
    for name, digest in json.loads(manifest.read_text())["sources"].items():
        original.extend(_classic_sections(str(source / name), digest))
    assert selected == [s for s in original if s["kind"] in ("elastic", "excitation")]
    assert len(selected) == 6
    assert sum(s["kind"] == "elastic" for s in selected) == 3
    assert [s["threshold_eV"] for s in selected if s["kind"] == "excitation"] == [
        2.39071, 3.5676, 1.177,
    ]


def test_realization_refuses_changed_source_before_publishing(tmp_path, sources):
    source, database = sources
    local = tmp_path / "source"
    local.mkdir()
    for name in ("N_LXCat.txt", "N_LXCat_extra.txt"):
        (local / name).write_bytes((source / name).read_bytes() + b"\n")
    destination = tmp_path / "realization"
    with pytest.raises(FingerprintMismatch):
        realize_sources(local, database / "input/eedf/atomic_nitrogen_channel_map.json", destination)
    assert not destination.exists()


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "threshold", "inverse", "extra_source"])
def test_realization_refuses_incomplete_frozen_policy(tmp_path, sources, mutation):
    source, database = sources
    frozen = json.loads((database / "input/eedf/atomic_nitrogen_channel_map.json").read_text())
    if mutation == "missing":
        frozen["channels"].pop()
    elif mutation == "duplicate":
        frozen["channels"].append(copy.deepcopy(frozen["channels"][0]))
    elif mutation == "threshold":
        frozen["channels"][1]["threshold_eV"] *= 1.01
    elif mutation == "inverse":
        del frozen["channels"][1]["inverse"]
    else:
        frozen["sources"]["unrelated.txt"] = "0" * 64
    path = tmp_path / "map.json"
    path.write_text(json.dumps(frozen))
    with pytest.raises(SpecError):
        realize_sources(source, path, tmp_path / "realization")


def test_refresh_preserves_existing_map_and_installs_complete_atomic_axes(tmp_path, sources):
    source, database = sources
    base = os.environ.get("EEDF_ATOMIC_BASE_SPEC")
    if not base:
        pytest.skip("realized Ar/N2 base spec is required")
    original = load_spec(base)
    original_channels = _load_channel_map(original["channel_map"])
    result = refresh_spec(base, database, source, tmp_path / "converted")
    spec = load_spec(result["spec_path"])
    channels = _load_channel_map(spec["channel_map"])
    assert channels[:-6] == original_channels
    assert spec["atomic_nitrogen_population"] == CONVENTION
    assert all(spec["axes"][axis] == [0.0, top] for axis, top in zip(AXES, BOUNDS))
    assert result["node_count"] == result["base_node_count"] * 8
    assert len(result["paired_owners"]) == 3
    assert {p[side] for p in result["paired_owners"] for side in ("ine", "sup")} == set(range(12, 18))
    assert spec["argon_4s_population"] == original["argon_4s_population"]
    assert spec["axes"]["f_Ar4s"] == original["axes"]["f_Ar4s"]
    assert result['zero_target_node_check']['passed']
    from rmgpy.tools.eedf.integrity import _resolved_properties
    point = {name: axis[0] for name, axis in spec['axes'].items() if name != 'u'}
    point.update({'f_Ar4s': 1e-5, **dict(zip(AXES, BOUNDS))})
    _, _, gas, states = _resolved_properties(spec, point)
    assert gas['N'] == pytest.approx(sum(BOUNDS))
    assert gas['Ar'] == pytest.approx(1 - point['x_N2'] - sum(BOUNDS))
    assert sum(states[name] for name in ('Ar(3P2)', 'Ar(3P1)', 'Ar(3P0)', 'Ar(1P1)')) / states['Ar(1S0)'] == pytest.approx(1e-5)
    pin = spec["input_files"]["Nitrogen/N-atomic.txt"]
    assert hashlib.sha256(Path(pin["path"]).read_bytes()).hexdigest() == pin["sha256"]
