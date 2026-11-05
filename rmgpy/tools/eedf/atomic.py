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

"""Realize pinned atomic nitrogen collisions and extend an Ar/N2 table spec.

Only the three elastic blocks and three reversible excitation blocks are
selected. Original headers and numerical rows are retained byte for byte;
atomic ionization remains an explicitly documented chemistry scope gap.
"""

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import re
import tempfile

from rmgpy.tools.eedf.atomic_population import (
    AXES, BOUNDS, CONVENTION, STATES, atomic_held_out_points,
)
from rmgpy.tools.eedf.channel_map import (
    _frozen_rows, _load_library, _owner_problem, _physical_channel_from_frozen,
    _read_pinned_bytes,
)
from rmgpy.tools.eedf.channels import (
    _classic_sections, _gas_values, _load_channel_map, _validate_physical_map,
)
from rmgpy.tools.eedf.integrity import _physical_properties
from rmgpy.tools.eedf.schema import FingerprintMismatch, SpecError, load_spec

COLLISION_SET = "IST-Lisbon-N-Wang2014"
LIBRARY = "PlasmaAtomicNitrogen"
SOURCE_NAMES = {"N_LXCat.txt", "N_LXCat_extra.txt"}
PAIRS = ((STATES[0], STATES[1], 2.39071),
         (STATES[0], STATES[2], 3.5676),
         (STATES[1], STATES[2], 1.177))
WEIGHTS = dict(zip(STATES, (4, 10, 6)))
EXPECTED = {
    **{f"e + {state} -> e + {state}, Elastic": 0.0 for state in STATES},
    **{f"e + {lower} <-> e + {upper}, Excitation": threshold
       for lower, upper, threshold in PAIRS},
}


def _read_frozen(path):
    raw = Path(path).read_bytes()
    frozen = json.loads(raw)
    rows = _frozen_rows(frozen, (LIBRARY,))
    if (frozen.get("set") != COLLISION_SET
            or set(frozen.get("sources", {})) != SOURCE_NAMES
            or set(frozen["loki_files"]) != {"N-atomic.txt"}
            or set(rows) != set(EXPECTED)):
        raise SpecError("atomic map requires the six sourced elastic/excitation blocks")
    for process, row in rows.items():
        excitation = " <-> " in process
        if (row["threshold_eV"] != EXPECTED[process]
                or row["class"] != ("A" if excitation else "B")
                or (excitation and row.get("inverse", {}).get("status") != "paired")):
            raise SpecError("atomic map threshold/class/paired inverse: " + process)
    return frozen, rows, hashlib.sha256(raw).hexdigest()


def realize_sources(source_directory, frozen_map, destination):
    """Write the frozen six-block selection only after verifying every source.

    The frozen database map pins both full downloaded files and the selected
    solver input. No numerical cross sections are redistributed in the repo.
    """
    frozen, rows, map_hash = _read_frozen(frozen_map)
    blocks, found = [], set()
    for name, digest in frozen["sources"].items():
        raw = _read_pinned_bytes(Path(source_directory) / name, digest).decode("utf-8")
        starts = list(re.finditer(r"(?m)^(ELASTIC|EXCITATION|IONIZATION)\s*$", raw))
        for index, start in enumerate(starts):
            end = starts[index + 1].start() if index + 1 < len(starts) else len(raw)
            block = raw[start.start():end]
            delimiters = list(re.finditer(r"(?m)^-{5,}\s*$", block))
            if len(delimiters) < 2:
                raise SpecError("incomplete atomic collision block")
            block = block[:delimiters[1].end()].strip() + "\n"
            match = re.search(r"COMMENT: \[(.*?)\]", block)
            if not match:
                raise SpecError("atomic collision lacks a process identity")
            process = match[1]
            if process not in rows:
                if start[1] != "IONIZATION":
                    raise SpecError("unmapped atomic collision: " + process)
                continue
            if process in found:
                raise SpecError("duplicate atomic collision: " + process)
            found.add(process)
            blocks.append(block)
    if found != set(EXPECTED):
        raise SpecError("atomic source selection is incomplete")
    content = "\n".join(blocks).encode("utf-8")
    digest = hashlib.sha256(content).hexdigest()
    if digest != frozen["loki_files"]["N-atomic.txt"]:
        raise FingerprintMismatch("atomic realized collision file")
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    path = destination / "N-atomic.txt"
    path.write_bytes(content)
    result = {"path": str(path), "sha256": digest, "kind": "cross_section"}
    (destination / "realization.json").write_text(json.dumps({
        "schema": "atomic-nitrogen-realization/1", "source_sha256": frozen["sources"],
        "frozen_map_sha256": map_hash, "realized": result,
        "selected_processes": list(rows), "excluded": "atomic ionization",
        "source_doi": "10.1103/PhysRevA.89.062714",
    }, indent=2) + "\n")
    return result


def refresh_spec(base_spec, database, source_directory, destination):
    """Preserve the base collision map and append three ordinary owner pairs.

    The existing Ar4s population policy, axes and reviewed inverse exceptions
    are retained. This creates a specification, never a production table.
    """
    # Reuse the reviewed population-owner validation after its upstream merge.
    from rmgpy.tools.eedf.channel_map import _paired_inverse

    base_spec, database = Path(base_spec).resolve(), Path(database).resolve()
    spec = load_spec(base_spec)
    base_count = math.prod(map(len, spec["axes"].values()))
    if set(spec["arm"]["gases"]) != {"Ar", "N2"} or set(AXES).intersection(spec["axes"]):
        raise SpecError("atomic refresh requires an Ar/N2 base without atomic axes")
    frozen_map = database / "input/eedf/atomic_nitrogen_channel_map.json"
    frozen, frozen_rows, map_hash = _read_frozen(frozen_map)
    channels = _load_channel_map(spec["channel_map"])
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    pin = realize_sources(source_directory, frozen_map, destination / "realization")
    spec["input_files"]["Nitrogen/N-atomic.txt"] = pin
    spec["arm"]["gases"].append("N")
    spec["gas_properties"]["fraction"].append("N = 0")
    spec["atomic_nitrogen_population"] = copy.deepcopy(CONVENTION)
    spec["axes"].update({name: [0.0, top] for name, top in zip(AXES, BOUNDS)})
    spec["refinement"]["max_rounds"] = 0
    spec["interpolation_accuracy"] = {
        name: {"rtol": 0.01, "atol": 0.0} for name in ("mobility_N", "field_power")
    }
    spec["screen"]["candidates"].extend(
        {"name": name, "reference": 0.0, "values": [0.0, top]}
        for name, top in zip(AXES, BOUNDS)
    )
    properties = spec["state_properties"]
    for key, values in (("population", (1.0, 0.0, 0.0)),
                        ("energy", (0.0, 2.39071, 3.5676)),
                        ("statisticalWeight", (4, 10, 6))):
        if not isinstance(properties.get(key), list):
            raise SpecError("atomic refresh requires inline state properties")
        if any(line.startswith("N(") for line in properties[key]):
            raise SpecError("atomic states already present in base properties")
        properties[key].extend(f"{state} = {value}" for state, value in zip(STATES, values))
    spec.setdefault("engine_state_map", []).extend(
        {"formula": "N", "electronic_state": term, "vibrational_level": -1,
         "multiplicity": multiplicity, "loki_state": state, "statistical_weight": weight}
        for state, term, multiplicity, weight in zip(STATES, ("", "2D", "2P"), (4, 2, 2), (4, 10, 6))
    )
    off_centre = [dict(point, **{name: 0.0 for name in AXES})
                  for point in spec["held_out"].get("off_centre", [])]
    off_centre.extend(atomic_held_out_points(spec["axes"]))
    spec["held_out"]["off_centre"] = off_centre
    gas, _ = _physical_properties(spec)
    mass_pin = spec["input_files"][gas["mass"]]
    masses = _gas_values(mass_pin["path"], mass_pin["sha256"])
    paired = []
    with tempfile.TemporaryDirectory(prefix="atomic-owner-") as directory:
        library, library_hashes = _load_library(database, Path(directory), LIBRARY)
        for physical in _classic_sections(pin["path"], pin["sha256"]):
            row = frozen_rows[physical["description"]]
            channel = _physical_channel_from_frozen(
                physical, row, COLLISION_SET, mass_values=masses, weights=WEIGHTS,
            )
            if row["class"] == "A":
                problem = _owner_problem(row, library)
                if problem:
                    raise SpecError("atomic ine owner: " + problem)
                pair = _paired_inverse(row, channel, library, COLLISION_SET)
                if pair is None:
                    raise SpecError("atomic excitation requires an ordinary inverse owner")
                paired.append(pair)
            channels.append(channel)
    _validate_physical_map(channels, spec, {name: axis[0] for name, axis in spec["axes"].items() if name != "u"})
    from rmgpy.tools.eedf.argon_population import check_zero_target_nodes
    population_nodes_checked = check_zero_target_nodes(spec, channels)
    channel_path = destination / "channel-map.json"
    content = (json.dumps(channels, indent=2) + "\n").encode()
    channel_path.write_bytes(content)
    spec["channel_map"] = {"path": str(channel_path), "sha256": hashlib.sha256(content).hexdigest()}
    spec_path = destination / "atomic-spec.json"
    # Validate in memory before publishing a runnable spec.
    from rmgpy.tools.eedf.schema import _validate_spec
    _validate_spec(spec)
    spec_path.write_text(json.dumps(spec, indent=2) + "\n")
    report = {"outcome": "converted", "spec_path": str(spec_path),
              "base_spec_sha256": hashlib.sha256(base_spec.read_bytes()).hexdigest(),
              "frozen_map_sha256": map_hash, "source_sha256": frozen["sources"],
              "library_sha256": library_hashes, "paired_owners": paired,
              "zero_target_node_check": {"passed": True, "population_nodes_checked": population_nodes_checked},
              "base_node_count": base_count, "node_count": base_count * 8, "cost_factor": 8,
              "scope": "spec only; no production table generated or qualified"}
    (destination / "conversion-report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-spec", required=True)
    parser.add_argument("--database", required=True)
    parser.add_argument("--sources", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    print(json.dumps(refresh_spec(args.base_spec, args.database, args.sources, args.output), indent=2))


if __name__ == "__main__":
    main()
