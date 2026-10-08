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

"""Acceptance tests for the canonical EEDF corroboration setup."""

import hashlib
import json
import math
import os
from pathlib import Path
import signal
import stat
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from unittest import mock

import h5py
import numpy as np
from scipy.constants import electron_mass

from rmgpy.tools.eedf.corroborate import (
    CorroborationError,
    ExecutionBackendUnavailable,
    _parse_bolos_blocks,
    _quadrature,
    build_canonical_setup,
    corroborate_row,
    SystemdScopeBackend,
)


class CanonicalSetupRedTest(unittest.TestCase):
    """The three r210 P1 bypasses must refuse at the production entry point."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.artifact = Path(self.temporary.name)
        self.collision = self.artifact / "argon.txt"
        self.mass = self.artifact / "masses.txt"
        self.mass.write_text("Ar %.17g\n" % (electron_mass / 1.0e-5))
        self._write_collision(1.0e-5)
        self.channel = {
            "description": "e + Ar(1S0) -> e + Ar(1S0), Elastic",
            "kind": "elastic",
            "classification": "B",
            "threshold_eV": 0.0,
            "target_fraction": 1.0,
            "product_fraction": 0.0,
            "sigma_max_m2": 1.0,
            "flux_group": "elastic",
            "reaction": None,
            "cross_section": {
                "energy_eV": [0.0, 1.0, 2.0],
                "sigma_m2": [1.0, 1.0, 1.0],
            },
            "mass_ratio": 1.0e-5,
        }
        self.channels = [self.channel]
        self.spec = {
            "arm": {
                "gases": ["Ar"],
                "additive": None,
                "feed_partial_pressure_Pa": 0.0,
            },
            "Tg_K": 298.15,
            "P_Pa": 100.0,
            "input_files": self._input_files(),
            "envelopes": {},
            "gas_properties": {"mass": "masses.txt", "fraction": ["Ar = 1.0"]},
            "state_properties": {"population": ["Ar(1S0) = 1.0"]},
            "working_conditions": {"excitationFrequency": 0},
            "solver_options": {
                "eedfType": "boltzmann",
                "ionizationOperatorType": "equalSharing",
                "growthModelType": "temporal",
                "includeEECollisions": False,
                "LXCatFiles": ["argon.txt"],
                "numerics": {
                    "energyGrid": {"maxEnergy": 2.0, "cellNumber": 2},
                    "maxPowerBalanceRelError": 1.0e-6,
                    "nonLinearRoutines": {"maxEedfRelError": 1.0e-7},
                },
            },
            "timeout_s": 1.0,
            "loki_commit": "a" * 40,
            "binary": {"path": "unused", "sha256": "b" * 64},
        }
        self.manifest = {
            "axes": {"u": [math.log(10.0)]},
            "energy_eV": [0.5, 1.5],
            "energy_edges_eV": [0.0, 1.0, 2.0],
            "branches": ["branch_0"],
            "channel_map": self.channels,
            "floors": {"absolute_power_share": 1.0e-4},
            "row_inputs": {"Tg_K": 298.15, "P_Pa": 100.0},
        }
        self._write_effective_config()
        self._write_metadata()
        self.f0 = np.exp(-np.array([0.5, 1.5]))
        self.f0 /= np.sum(np.sqrt([0.5, 1.5]) * self.f0)
        self._write_table(np.array([1.0]))
        self.marker = self.artifact / "solver-ran"
        self.solver_path = self.artifact / "fake-solver.py"
        self._write_solver()

    def tearDown(self):
        self.temporary.cleanup()

    def _write_collision(self, header_mass_ratio):
        self.collision.write_text(
            "ELASTIC\n"
            "Ar\n"
            f"{header_mass_ratio:.17g}\n"
            "PROCESS: E + Ar -> E + Ar, Elastic\n"
            f"PARAM.: m/M = {header_mass_ratio:.17g}\n"
            "[e + Ar(1S0) -> e + Ar(1S0), Elastic]\n"
            "----------------\n"
            "0 1\n1 1\n2 1\n"
            "----------------\n"
        )

    def _input_files(self):
        return {
            "argon.txt": {
                "path": str(self.collision),
                "sha256": hashlib.sha256(self.collision.read_bytes()).hexdigest(),
                "kind": "cross_section",
            },
            "masses.txt": {
                "path": str(self.mass),
                "sha256": hashlib.sha256(self.mass.read_bytes()).hexdigest(),
                "kind": "property",
            },
        }

    def _write_metadata(self):
        (self.artifact / "generation_spec.json").write_text(json.dumps(self.spec))
        (self.artifact / "manifest.json").write_text(json.dumps(self.manifest))

    def _write_effective_config(self):
        path = self.artifact / "loki-effective.json"
        path.write_text(
            json.dumps(
                {
                    "nativeGasProperties": {
                        "fraction": {"Ar": 1.0},
                        "mass": {"Ar": electron_mass / 1.0e-5, "e": electron_mass},
                    },
                    "nativeCollisionDescriptions": [
                        channel["description"] for channel in self.channels
                    ],
                    "workingConditions": {
                        "gasTemperature": {"unit": "K", "value": 298.15},
                        "gasPressure": {"unit": "Pa", "value": 100.0},
                        "reducedField": {"unit": "Td", "value": 10.0},
                    },
                    "electronKinetics": {
                        "gasProperties": {"fraction": {"Ar": 1.0}},
                        "stateProperties": {
                            "population": {
                                "states": {"Ar(1S0)": {"type": "constant", "value": 1.0}}
                            }
                        },
                        "ionizationOperatorType": "equalSharing",
                        "growthModelType": "temporal",
                        "includeEECollisions": False,
                        "LXCatFiles": ["argon.txt"],
                        "numerics": {
                            "energyGrid": {"maxEnergy": 2.0, "cellNumber": 2},
                            "nonLinearRoutines": {"maxEedfRelError": 1.0e-7},
                        },
                    },
                }
            )
        )
        self.manifest["loki_effective_config"] = {
            "path": path.name,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }

    def _write_table(self, target_fractions):
        row = {
            "energy_eV": np.array([0.5, 1.5]),
            "energy_edges_eV": np.array([0.0, 1.0, 2.0]),
            "gas_temperature_K": 298.15,
            "target_fractions": target_fractions,
            "product_fractions": np.zeros(len(self.channels)),
        }
        processes = [
            {
                **channel,
                "threshold": channel["threshold_eV"],
                "energy": channel["cross_section"]["energy_eV"],
                "sigma": channel["cross_section"]["sigma_m2"],
            }
            for channel in self.channels
        ]
        values = _quadrature(row["energy_eV"], self.f0, processes, row)
        with h5py.File(self.artifact / "table.h5", "w") as h5:
            group = h5.create_group("branches").create_group("branch_0")
            group.create_dataset("f0", data=[self.f0])
            group.create_dataset("converged", data=[True])
            group.create_dataset(
                "termination_status", data=["converged"], dtype=h5py.string_dtype()
            )
            group.create_dataset("convergence_residual", data=[1.0e-9])
            group.create_dataset("iteration_count", data=[1])
            group.create_dataset("convergence_tolerance", data=[1.0e-7])
            group.create_dataset(
                "channels",
                data=[[channel["description"] for channel in self.channels]],
                dtype=h5py.string_dtype(),
            )
            group.create_dataset("target_fractions", data=[target_fractions])
            group.create_dataset("product_fractions", data=[np.zeros(len(self.channels))])
            group.create_dataset("rate_floors", data=[[1.0e-40] * len(self.channels)])
            group.create_dataset(
                "below_floor", data=[[[False, False] for _ in self.channels]]
            )
            group.create_dataset("k_ine", data=[values["rates"]])
            group.create_dataset("channel_power", data=[values["channel_power"]])
            swarm = group.create_group("swarm")
            swarm.create_dataset("mean_energy_eV", data=[values["mean_energy_eV"]])
            swarm.create_dataset("mobility_N", data=[values["mobility_N"]])
            power = group.create_group("power_groups")
            power.create_dataset("field", data=[max(abs(values["total_power"]), 1.0)])
        evidence = self.artifact / "loki-convergence.json"
        native_log = self.artifact / "loki-convergence.log"
        native_log.write_text(
            "T9_CONVERGENCE mode=linear converged=1 iterations=1 "
            "residual=1e-9 tolerance=1e-7\n"
        )
        evidence.write_text(
            json.dumps(
                {
                    "mode": "linear",
                    "termination_status": "converged",
                    "converged": True,
                    "residual": 1.0e-9,
                    "iteration_count": 1,
                    "tolerance": 1.0e-7,
                    "provenance": {
                        "base_solver_commit": "a" * 40,
                        "base_solver_sha256": "b" * 64,
                        "evidence_log_sha256": hashlib.sha256(native_log.read_bytes()).hexdigest(),
                        "instrumentation_source_sha256": "c" * 64,
                        "instrumented_solver_sha256": "d" * 64,
                    },
                }
            )
        )
        self.manifest["loki_convergence_evidence"] = {
            "path": evidence.name,
            "sha256": hashlib.sha256(evidence.read_bytes()).hexdigest(),
        }
        self.manifest["loki_convergence_log"] = {
            "path": native_log.name,
            "sha256": hashlib.sha256(native_log.read_bytes()).hexdigest(),
        }
        self._write_metadata()

    def _write_solver(self):
        result = {
            "status": "converged",
            "energy": [0.5, 1.5],
            "f0": self.f0.tolist(),
            "version": "fake",
        }
        self.solver_path.write_text(
            "#!/usr/bin/env python3\n"
            "from pathlib import Path\n"
            f"Path({str(self.marker)!r}).write_text('ran')\n"
            f"print({json.dumps(result)!r})\n"
        )
        self.solver_path.chmod(self.solver_path.stat().st_mode | stat.S_IXUSR)

    def _add_ionization(self):
        with self.collision.open("a") as stream:
            stream.write(
                "IONIZATION\n"
                "Ar -> Ar+\n"
                "0.5\n"
                "PROCESS: E + Ar -> E + E + Ar+, Ionization\n"
                "PARAM.: E = 0.5 eV\n"
                "[e + Ar(1S0) -> e + e + Ar(+,gnd), Ionization]\n"
                "----------------\n"
                "0 0\n1 0.25\n2 0.25\n"
                "----------------\n"
            )
        self.channels.append(
            {
                "description": "e + Ar(1S0) -> e + e + Ar(+,gnd), Ionization",
                "kind": "ionization",
                "classification": "B",
                "threshold_eV": 0.5,
                "target_fraction": 1.0,
                "product_fraction": 0.0,
                "sigma_max_m2": 0.25,
                "flux_group": "ionization",
                "reaction": None,
                "cross_section": {
                    "energy_eV": [0.0, 1.0, 2.0],
                    "sigma_m2": [0.0, 0.25, 0.25],
                },
                "opb_eV": 0.5,
            }
        )
        self.spec["input_files"] = self._input_files()
        self.manifest["channel_map"] = self.channels
        self._write_effective_config()
        self._write_metadata()
        self._write_table(np.ones(len(self.channels)))

    def _bolos_result(self, *, overrides=None, omit=()):
        ionizing = len(self.channels) == 2
        processes = [
            {
                "kind": "elastic",
                "target": "Ar",
                "product": None,
                "threshold": 0.0,
                "mass_ratio": 1.0e-5,
                "energy": [0.0, 1.0, 2.0],
                "sigma": [1.0, 1.0, 1.0],
            }
        ]
        if ionizing:
            processes.append(
                {
                    "kind": "ionization",
                    "target": "Ar",
                    "product": "Ar+",
                    "threshold": 0.5,
                    "mass_ratio": None,
                    "energy": [0.0, 1.0, 2.0],
                    "sigma": [0.0, 0.25, 0.25],
                }
            )
        result = {
            "termination_status": "converged",
            "converged": True,
            "iterations": 7,
            "residual": 1.0e-8,
            "tolerance": 1.0e-7,
            "energy": [0.5, 1.5],
            "f0": self.f0.tolist(),
            "native_mobility_N": 101136.53988549036 if ionizing else 143932.90827901874,
            "native_mean_energy_eV": 0.8891958083193515,
            "version": "fake-native",
            "effective_config": {
                "Tg_K": 298.15,
                "EN_Td": 10.0,
                "grid": {"max_energy_eV": 2.0, "cell_count": 2},
                "gas_fractions": {"Ar": 1.0},
                "target_mass_ratios": {"Ar": 1.0e-5},
                "processes": processes,
                "ionization_energy_sharing": "equalSharing",
                "electron_growth": "temporal",
                "mobility_definition": "temporal-growth-corrected reduced mobility",
            },
        }
        result.update(overrides or {})
        for name in omit:
            result.pop(name, None)
        return result

    def _execution(self, result=None, **overrides):
        execution = {
            "returncode": 0,
            "stdout": json.dumps(result or self._bolos_result()) + "\n",
            "stderr": "",
            "cleanup_passed": True,
            "backend": "systemd-user-scope-cgroup.kill",
        }
        execution.update(overrides)
        return execution

    def _run(self, node=None):
        return corroborate_row(
            self.artifact,
            "branch_0",
            node or {},
            math.log(10.0),
            solvers={
                "bolos": {
                    "kind": "bolos",
                    "interpreter": str(self.solver_path),
                }
            },
        )

    def test_bolos_header_mass_mismatch_refuses_before_solver_runs(self):
        self._write_collision(2.0e-5)
        self.spec["input_files"] = self._input_files()
        self._write_metadata()
        with self.assertRaisesRegex(CorroborationError, "mass.*mismatch"):
            self._run()
        self.assertFalse(self.marker.exists())

    def test_bolos_header_participant_must_match_loki_identity(self):
        self._write_collision(1.0e-5)
        self.collision.write_text(self.collision.read_text().replace("Ar\n1.000", "Xe\n1.000", 1))
        self.spec["input_files"] = self._input_files()
        self._write_metadata()
        with self.assertRaisesRegex(CorroborationError, "participant mismatch"):
            self._run()
        self.assertFalse(self.marker.exists())

    def test_caller_temperature_override_refuses_before_solver_runs(self):
        with self.assertRaisesRegex(CorroborationError, "Tg_K.*override|override.*Tg_K"):
            self._run({"Tg_K": 301.0})
        self.assertFalse(self.marker.exists())

    def test_stored_fractions_must_match_pinned_populations_before_solver_runs(self):
        self._write_table(np.array([0.5]))
        with self.assertRaisesRegex(CorroborationError, "fraction.*population|population.*fraction"):
            self._run()
        self.assertFalse(self.marker.exists())

    def test_hand_audited_parser_fixture_has_expected_physical_values(self):
        blocks = _parse_bolos_blocks(self.collision.read_bytes(), self.collision.name)
        self.assertEqual(
            blocks,
            [
                {
                    "kind": "elastic",
                    "target": "Ar",
                    "product": None,
                    "threshold": 0.0,
                    "mass_ratio": 1.0e-5,
                    "energy": (0.0, 1.0, 2.0),
                    "sigma": (1.0, 1.0, 1.0),
                    "comments": (
                        "PROCESS: E + Ar -> E + Ar, Elastic",
                        "PARAM.: m/M = 1.0000000000000001e-05",
                        "[e + Ar(1S0) -> e + Ar(1S0), Elastic]",
                    ),
                }
            ],
        )

    def test_canonical_setup_is_immutable_and_version_fingerprinted(self):
        setup = build_canonical_setup(
            self.artifact, "branch_0", {}, math.log(10.0)
        )
        self.assertEqual(setup.fingerprint_version, "eedf-canonical-setup-v1")
        self.assertRegex(setup.fingerprint, r"^[0-9a-f]{64}$")
        with self.assertRaises(AttributeError):
            setup.gas_temperature_K = 301.0

    def test_changed_pinned_source_refuses_before_solver_runs(self):
        self.collision.write_text(self.collision.read_text() + "\n")
        with self.assertRaisesRegex(CorroborationError, "SHA-256 mismatch"):
            self._run()
        self.assertFalse(self.marker.exists())

    def test_loki_native_effective_fraction_mismatch_refuses_before_solver_runs(self):
        path = self.artifact / "loki-effective.json"
        effective = json.loads(path.read_text())
        effective["nativeGasProperties"]["fraction"]["Ar"] = 0.5
        path.write_text(json.dumps(effective))
        self.manifest["loki_effective_config"]["sha256"] = hashlib.sha256(
            path.read_bytes()
        ).hexdigest()
        self._write_metadata()
        with self.assertRaisesRegex(CorroborationError, "gas_fractions"):
            self._run()
        self.assertFalse(self.marker.exists())

    def test_loki_native_units_and_population_type_are_semantic(self):
        for field, value, diagnostic in (
            ("unit", "eV", "Tg_K"),
            ("type", "dynamic", "state_populations"),
        ):
            with self.subTest(field=field):
                self._write_effective_config()
                path = self.artifact / "loki-effective.json"
                effective = json.loads(path.read_text())
                if field == "unit":
                    effective["workingConditions"]["gasTemperature"][field] = value
                else:
                    effective["electronKinetics"]["stateProperties"]["population"]["states"]["Ar(1S0)"][field] = value
                path.write_text(json.dumps(effective))
                self.manifest["loki_effective_config"]["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
                self._write_metadata()
                with self.assertRaisesRegex(CorroborationError, diagnostic):
                    self._run()
                self.assertFalse(self.marker.exists())

    def test_loki_false_convergence_refuses(self):
        with h5py.File(self.artifact / "table.h5", "r+") as h5:
            h5["branches/branch_0/converged"][0] = False
        with self.assertRaisesRegex(CorroborationError, "converged evidence is false"):
            self._run()

    def test_loki_missing_convergence_evidence_refuses(self):
        with h5py.File(self.artifact / "table.h5", "r+") as h5:
            del h5["branches/branch_0/iteration_count"]
        with self.assertRaisesRegex(CorroborationError, "evidence is missing: iteration_count"):
            self._run()

    def test_loki_malformed_convergence_evidence_refuses(self):
        with h5py.File(self.artifact / "table.h5", "r+") as h5:
            group = h5["branches/branch_0"]
            del group["converged"]
            group.create_dataset("converged", data=["False"], dtype=h5py.string_dtype())
        with self.assertRaisesRegex(CorroborationError, "not boolean"):
            self._run()

    def test_loki_nonfinite_convergence_residual_refuses(self):
        with h5py.File(self.artifact / "table.h5", "r+") as h5:
            h5["branches/branch_0/convergence_residual"][0] = np.nan
        with self.assertRaisesRegex(CorroborationError, "residual.*finite"):
            self._run()

    def test_numeric_strings_are_not_silently_coerced(self):
        with h5py.File(self.artifact / "table.h5", "r+") as h5:
            group = h5["branches/branch_0"]
            del group["target_fractions"]
            group.create_dataset("target_fractions", data=[["1.0"]], dtype=h5py.string_dtype())
        with self.assertRaisesRegex(CorroborationError, "numeric array"):
            self._run()

    def test_loki_native_convergence_mismatch_refuses(self):
        path = self.artifact / "loki-convergence.json"
        evidence = json.loads(path.read_text())
        evidence["iteration_count"] = 2
        path.write_text(json.dumps(evidence))
        self.manifest["loki_convergence_evidence"]["sha256"] = hashlib.sha256(
            path.read_bytes()
        ).hexdigest()
        self._write_metadata()
        with self.assertRaisesRegex(CorroborationError, "differs from native log"):
            self._run()

    @mock.patch("rmgpy.tools.eedf.corroborate.SystemdScopeBackend.run")
    def test_bolos_false_missing_malformed_and_nonfinite_convergence_refuse(self, run):
        cases = (
            self._bolos_result(overrides={"converged": False}),
            self._bolos_result(omit=("residual",)),
            self._bolos_result(overrides={"converged": "True"}),
            self._bolos_result(overrides={"residual": float("nan")}),
        )
        for result in cases:
            with self.subTest(result=result):
                run.return_value = self._execution(result)
                with self.assertRaisesRegex(CorroborationError, "bolos"):
                    self._run()

    @mock.patch("rmgpy.tools.eedf.corroborate.SystemdScopeBackend.run")
    def test_bolos_effective_header_mass_mismatch_refuses(self, run):
        result = self._bolos_result()
        result["effective_config"]["target_mass_ratios"]["Ar"] = 2.0e-5
        run.return_value = self._execution(result)
        with self.assertRaisesRegex(CorroborationError, "parsed mass ratio"):
            self._run()

    @mock.patch("rmgpy.tools.eedf.corroborate.SystemdScopeBackend.run")
    def test_snapshot_is_used_when_source_changes_after_setup_creation(self, run):
        original = self.collision.read_bytes()

        def mutate_and_finish(argv, *, cwd, timeout_s):
            request = json.loads(Path(argv[-1]).read_text())
            self.collision.write_text("changed after canonical snapshot")
            self.assertEqual(Path(request["collision_files"][0]).read_bytes(), original)
            return self._execution()

        run.side_effect = mutate_and_finish
        result = self._run()
        self.assertEqual(result["verdict"], "PASS")

    @mock.patch("rmgpy.tools.eedf.corroborate.SystemdScopeBackend.run")
    def test_ionization_free_positive_shows_every_pass_term(self, run):
        run.return_value = self._execution()
        result = self._run()
        self.assertTrue(result["accepted"])
        self.assertEqual(result["verdict"], "PASS")
        self.assertTrue(all(result["pass_terms"].values()))
        self.assertEqual(result["effective_configuration_mismatches"], {"loki": [], "bolos": []})

    @mock.patch("rmgpy.tools.eedf.corroborate.SystemdScopeBackend.run")
    def test_ionizing_positive_has_nonzero_growth_contribution(self, run):
        self._add_ionization()
        run.return_value = self._execution()
        result = self._run()
        self.assertEqual(result["verdict"], "PASS")
        self.assertGreater(result["transport_definition"]["growth_rate_m3_s"], 0.0)

    @mock.patch("rmgpy.tools.eedf.corroborate.SystemdScopeBackend.run")
    def test_ionizing_comparison_omitting_growth_fails_transport_definition(self, run):
        self._add_ionization()
        result = self._bolos_result(overrides={"native_mobility_N": 115146.32662321498})
        run.return_value = self._execution(result)
        with self.assertRaisesRegex(CorroborationError, "transport-definition mismatch"):
            self._run()

    @mock.patch("rmgpy.tools.eedf.corroborate.SystemdScopeBackend.run")
    def test_surviving_child_cleanup_failure_emits_no_result(self, run):
        run.side_effect = CorroborationError(
            "solver execution cleanup failed; child survived"
        )
        with self.assertRaisesRegex(CorroborationError, "child survived"):
            self._run()

    def test_caller_cannot_loosen_frozen_tolerance(self):
        with self.assertRaisesRegex(CorroborationError, "tolerance is frozen"):
            corroborate_row(
                self.artifact,
                "branch_0",
                {},
                math.log(10.0),
                solvers={"bolos": {"kind": "bolos", "interpreter": "unused"}},
                rtol=0.1,
            )

    @mock.patch("rmgpy.tools.eedf.corroborate.shutil.which", return_value=None)
    def test_unavailable_descendant_backend_refuses_with_named_reason(self, which):
        with self.assertRaisesRegex(ExecutionBackendUnavailable, "systemd user-scope"):
            SystemdScopeBackend()


@unittest.skipUnless(
    os.environ.get("EEDF_RUN_CGROUP_TESTS") == "1",
    "requires the deployment host's delegated user systemd scope",
)
class SystemdScopeLifecycleTest(unittest.TestCase):
    """Exercise descendant ownership and cgroup.kill on the deployment host."""

    def setUp(self):
        self.backend = SystemdScopeBackend()
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.unrelated = subprocess.Popen(["/bin/sleep", "30"])

    def tearDown(self):
        self.assertIsNone(self.unrelated.poll(), "backend touched an unrelated process")
        self.unrelated.terminate()
        self.unrelated.wait(timeout=2)
        self.temporary.cleanup()

    def test_normal_exit_leaves_no_job_process(self):
        result = self.backend.run(["/bin/sh", "-c", "exit 0"], cwd=self.root, timeout_s=2)
        self.assertEqual(result["returncode"], 0)
        self.assertTrue(result["cleanup_passed"])

    def test_solver_failure_leaves_no_job_process(self):
        result = self.backend.run(["/bin/sh", "-c", "exit 7"], cwd=self.root, timeout_s=2)
        self.assertEqual(result["returncode"], 7)
        self.assertTrue(result["cleanup_passed"])

    def test_timeout_kills_job_process(self):
        marker = self.root / "timeout-survivor"
        command = "sleep 1; touch " + str(marker)
        with self.assertRaisesRegex(CorroborationError, "timed out; cleanup_passed=True"):
            self.backend.run(["/bin/sh", "-c", command], cwd=self.root, timeout_s=0.1)
        time.sleep(1.1)
        self.assertFalse(marker.exists())

    def test_interruption_kills_job_process(self):
        marker = self.root / "interrupt-survivor"
        command = "sleep 1; touch " + str(marker)
        timer = threading.Timer(0.2, os.kill, args=(os.getpid(), signal.SIGINT))
        timer.start()
        try:
            with self.assertRaises(KeyboardInterrupt):
                self.backend.run(["/bin/sh", "-c", command], cwd=self.root, timeout_s=2)
        finally:
            timer.cancel()
        time.sleep(1.1)
        self.assertFalse(marker.exists())

    def test_double_forked_descendant_is_killed(self):
        marker = self.root / "double-fork-survivor"
        helper = self.root / "double_fork.py"
        helper.write_text(
            "import os, pathlib, time\n"
            "if os.fork() == 0:\n"
            "    os.setsid()\n"
            "    if os.fork() == 0:\n"
            "        time.sleep(1)\n"
            f"        pathlib.Path({str(marker)!r}).write_text('survived')\n"
            "        os._exit(0)\n"
            "    os._exit(0)\n"
        )
        with self.assertRaisesRegex(CorroborationError, "timed out; cleanup_passed=True"):
            self.backend.run([sys.executable, str(helper)], cwd=self.root, timeout_s=0.1)
        time.sleep(1.1)
        self.assertFalse(marker.exists())


if __name__ == "__main__":
    unittest.main()
