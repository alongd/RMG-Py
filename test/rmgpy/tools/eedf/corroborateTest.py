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
import io
import json
import math
import os
from pathlib import Path
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import threading
import time
import types
import unittest
from unittest import mock

import h5py
import numpy as np
import yaml
from scipy.constants import electron_mass

from rmgpy.tools.eedf import corroborate
from rmgpy.tools.eedf.corroborate import (
    _BOLOS_WORKER,
    _OWNER_WATCH,
    CorroborationError,
    ExecutionBackendUnavailable,
    _bolos_request,
    _loki_effective_configuration,
    _loki_mismatches_live,
    _parse_bolos_blocks,
    _quadrature,
    _same_effective_value,
    build_canonical_setup,
    corroborate_row,
    SystemdScopeBackend,
)


REPOSITORY = Path(__file__).resolve().parents[4]


def repository_env(**extra):
    """Environment for a Python subprocess that must import this checkout's rmgpy.

    The subprocess runs elsewhere (a temporary directory), so neither the
    caller's PYTHONPATH nor its working directory can be relied on to put this
    checkout ahead of another installed rmgpy.
    """
    paths = [str(REPOSITORY)] + [p for p in os.environ.get("PYTHONPATH", "").split(os.pathsep) if p]
    return dict(os.environ, PYTHONPATH=os.pathsep.join(paths), **extra)


def own_shutdown_signals(test):
    """Give SIGINT/SIGTERM Python's own dispositions for one test, then restore.

    Code that ran earlier in the process may have left a foreign C handler:
    RDKit's EmbedMolecule (reached by the export census through rd_embed)
    leaves a SIGINT handler that swallows the first SIGINT and then resets it
    to SIG_DFL. run() correctly leaves such a disposition alone, so a test that
    sends itself SIGINT and expects KeyboardInterrupt has to establish
    Python's handler first. The prior kernel dispositions are restored after.
    """
    saved = {}
    for signum, handler in ((signal.SIGINT, signal.default_int_handler), (signal.SIGTERM, signal.SIG_DFL)):
        saved[signum] = (signal.getsignal(signum), corroborate._raw_sigaction(signum))
        signal.signal(signum, handler)

    def restore():
        for signum, (view, raw) in saved.items():
            if view is not None:
                signal.signal(signum, view)
            if raw is not None:
                corroborate._restore_raw_sigaction(signum, raw)

    test.addCleanup(restore)


class CanonicalSetupRedTest(unittest.TestCase):
    """The three r210 P1 bypasses must refuse at the production entry point."""

    def setUp(self):
        own_shutdown_signals(self)
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
            "omp_threads": 1,
            "tolerances": {
                "normalization_atol": 1.0e-8,
                "fingerprint_rtol": 1.0e-12,
            },
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
        self.spec["binary"] = {
            "path": str(self.solver_path),
            "sha256": hashlib.sha256(self.solver_path.read_bytes()).hexdigest(),
        }
        self._write_metadata()
        self.loki_patcher = mock.patch(
            "rmgpy.tools.eedf.corroborate._run_loki"
        )
        self.run_loki = self.loki_patcher.start()
        self.run_loki.side_effect = lambda setup, config: (
            self._loki_result(), self._execution(self._loki_result())
        )
        self.addCleanup(self.loki_patcher.stop)

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
                "energy_intervals": [[1.0e-9, 1.0], [1.0, 1.999999999]],
                "sigma_intervals": [[1.0, 1.0], [1.0, 1.0]],
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
                    "energy_intervals": [[0.5, 1.0], [1.0, 1.999999999]],
                    "sigma_intervals": [[0.125, 0.25], [0.25, 0.25]],
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
                "populations": {"Ar": 1.0},
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

    def _loki_result(self, *, overrides=None):
        result = self._bolos_result()
        result["termination_status"] = "converged"
        result["converged"] = True
        result["residual"] = 1.0e-8
        result["tolerance"] = 1.0e-6
        result["version"] = "fake-loki"
        result.update(overrides or {})
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
                "loki": {
                    "kind": "loki",
                    "executable": str(self.solver_path),
                },
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
        setup = build_canonical_setup(
            self.artifact, "branch_0", {}, math.log(10.0)
        )
        self.assertEqual(setup.target_fractions, (1.0,))

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

    def test_stored_loki_sidecars_are_not_live_adapter_inputs(self):
        for name in ("loki-effective.json", "loki-convergence.json", "loki-convergence.log"):
            (self.artifact / name).unlink()
        for name in ("loki_effective_config", "loki_convergence_evidence", "loki_convergence_log"):
            self.manifest.pop(name)
        self._write_metadata()
        setup = build_canonical_setup(
            self.artifact, "branch_0", {}, math.log(10.0)
        )
        self.assertRegex(setup.fingerprint, r"^[0-9a-f]{64}$")

    def test_bolos_request_uses_canonical_state_population(self):
        self.spec["state_properties"]["population"] = ["Ar(1S0) = 0.999"]
        self.channel["target_fraction"] = 0.999
        self._write_metadata()
        self._write_table(np.array([0.999]))
        effective_path = self.artifact / "loki-effective.json"
        effective = json.loads(effective_path.read_text())
        effective["electronKinetics"]["stateProperties"]["population"]["states"]["Ar(1S0)"]["value"] = 0.999
        effective_path.write_text(json.dumps(effective))
        self.manifest["loki_effective_config"]["sha256"] = hashlib.sha256(
            effective_path.read_bytes()
        ).hexdigest()
        self._write_metadata()
        setup = build_canonical_setup(
            self.artifact, "branch_0", {}, math.log(10.0)
        )
        with tempfile.TemporaryDirectory() as temporary:
            request = json.loads(
                _bolos_request(setup, Path(temporary)).read_text()
            )
        self.assertEqual(request["populations"], {"Ar": 0.999})

    def test_effective_float_comparison_has_no_absolute_zero_floor(self):
        self.assertFalse(_same_effective_value(0.0, 1.0e-20))

    def test_live_loki_effective_config_covers_required_numerical_fields(self):
        setup = build_canonical_setup(
            self.artifact, "branch_0", {}, math.log(10.0)
        )
        with tempfile.TemporaryDirectory() as temporary:
            native_setup = Path(temporary) / "output" / "solve" / "setup.json"
            native_setup.parent.mkdir(parents=True)
            native = yaml.safe_load(setup.loki_setup)
            working = native["workingConditions"]
            for field, unit in (
                ("gasTemperature", "K"),
                ("gasPressure", "Pa"),
                ("excitationFrequency", "Hz"),
            ):
                working[field] = {"unit": unit, "value": working[field]}
            kinetics = native["electronKinetics"]
            kinetics["gasProperties"]["fraction"] = {"Ar": 1.0}
            kinetics["stateProperties"]["population"] = {
                "states": {"Ar(1S0)": {"type": "constant", "value": 1.0}}
            }
            native_setup.write_text(json.dumps(native))
            effective = _loki_effective_configuration(native_setup)
        self.assertEqual(_loki_mismatches_live(setup, effective), [])
        mutations = {
            "excitation_frequency_Hz": 1.0,
            "eedf_mode": "prescribed",
            "electron_electron_collisions": True,
            "convergence_tolerance": 2.0e-7,
        }
        for field, value in mutations.items():
            with self.subTest(field=field):
                changed = dict(effective)
                changed[field] = value
                self.assertTrue(
                    any(field in mismatch for mismatch in _loki_mismatches_live(setup, changed))
                )

    def test_bolos_worker_captures_initialized_process_arrays(self):
        target = types.SimpleNamespace(name="Ar", density=0.0, mass_ratio=1.0e-5)
        process = types.SimpleNamespace(
            kind="ELASTIC",
            product=None,
            threshold=0.0,
            mass_ratio=1.0e-5,
            data=np.array([[7.0, 70.0], [8.0, 80.0]]),
            eps=np.array([[0.0, 1.0], [1.0, 2.0]]),
            sigma=np.array([[3.0, 4.0], [4.0, 5.0]]),
        )

        class FakeSolver:
            def __init__(self, collision_grid):
                self.target = {"Ar": target}
                self.kT = self.EN = None
                self.benergy = np.array([0.0, 1.0, 2.0])
                self.cenergy = np.array([0.5, 1.5])
                self.sigma_m = np.ones(3)
                self.n = 2
                self.initialized = False

            def load_collisions(self, parsed):
                pass

            def init(self):
                self.initialized = True

            def iter_all(self):
                self.assert_initialized()
                return iter(((target, process),))

            def iter_growth(self):
                return iter((object(),))

            def assert_initialized(self):
                if not self.initialized:
                    raise AssertionError("worker read processes before init")

            def maxwell(self, mean_energy):
                return np.ones(2)

            def converge(self, initial, **kwargs):
                return np.array([0.5, 0.5]), 1, 1.0e-9

            def mobility(self, distribution):
                return 2.0

            def mean_energy(self, distribution):
                return 3.0

        bolos = types.ModuleType("bolos")
        bolos.__version__ = "test"
        bolos.grid = types.SimpleNamespace(LinearGrid=lambda *args: object())
        bolos.parser = types.SimpleNamespace(parse=lambda stream: [])
        bolos.solver = types.SimpleNamespace(
            BoltzmannSolver=FakeSolver,
            ConvergenceError=RuntimeError,
            TOWNSEND=1.0,
            GAMMA=1.0,
        )
        request = {
            "collision_files": [str(self.collision)],
            "populations": {"Ar": 0.999},
            "Tg_K": 300.0,
            "EN_Td": 10.0,
            "grid": {"max_energy_eV": 2.0, "cell_count": 2},
            "convergence": {"max_iterations": 10, "tolerance": 1.0e-7},
        }
        request_path = self.artifact / "worker-request.json"
        request_path.write_text(json.dumps(request))
        output = io.StringIO()
        with mock.patch.dict(sys.modules, {"bolos": bolos}), mock.patch.object(
            sys, "argv", ["bolos-worker", str(request_path)]
        ), mock.patch("sys.stdout", output):
            exec(_BOLOS_WORKER, {})
        emitted = json.loads(output.getvalue().splitlines()[-1])
        effective = emitted["effective_config"]["processes"][0]
        self.assertEqual(effective["energy_intervals"], process.eps.tolist())
        self.assertEqual(effective["sigma_intervals"], process.sigma.tolist())
        self.assertNotEqual(effective["energy_intervals"], process.data.tolist())

    def test_changed_pinned_source_refuses_before_solver_runs(self):
        self.collision.write_text(self.collision.read_text() + "\n")
        with self.assertRaisesRegex(CorroborationError, "SHA-256 mismatch"):
            self._run()
        self.assertFalse(self.marker.exists())

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
        self.run_loki.assert_called_once()
        run.assert_called_once()

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
                solvers={
                    "loki": {"kind": "loki", "executable": "unused"},
                    "bolos": {"kind": "bolos", "interpreter": "unused"},
                },
                rtol=0.1,
            )

    @mock.patch("rmgpy.tools.eedf.corroborate.shutil.which", return_value=None)
    def test_unavailable_descendant_backend_refuses_with_named_reason(self, which):
        with self.assertRaisesRegex(ExecutionBackendUnavailable, "systemd user-scope"):
            SystemdScopeBackend()

    def test_unknown_process_group_is_not_clean(self):
        self.assertFalse(SystemdScopeBackend._empty(None))

    def test_scope_runs_payload_behind_owner_watch(self):
        backend = object.__new__(SystemdScopeBackend)
        backend.systemd_run = "/usr/bin/systemd-run"
        backend.systemctl = "/usr/bin/systemctl"
        popen = mock.Mock(side_effect=OSError("no launcher"))
        with mock.patch("rmgpy.tools.eedf.corroborate.subprocess.Popen", popen):
            with self.assertRaisesRegex(ExecutionBackendUnavailable, "could not start"):
                backend.run(["solver", "--flag"], cwd=self.artifact, timeout_s=1.0)
        command = popen.call_args.args[0]
        unit = next(arg.split("=", 1)[1] for arg in command if arg.startswith("--unit="))
        watch = command.index(_OWNER_WATCH)
        lifeline = int(command[watch + 2])
        self.assertEqual(command[:4], ["/usr/bin/systemd-run", "--user", "--scope", "--quiet"])
        self.assertEqual(command[watch - 4:watch], [sys.executable, "-I", "-S", "-c"])
        self.assertEqual(
            command[watch + 1:],
            [str(os.getpid()), str(lifeline), unit, "/usr/bin/systemctl", "--", "solver", "--flag"],
        )
        # Only the lifeline's read end is passed; both ends are closed once run() is over.
        self.assertEqual(popen.call_args.kwargs["pass_fds"], (lifeline,))
        self.assertEqual(corroborate._LIFELINES, {})
        with self.assertRaises(OSError):
            os.fstat(lifeline)

    def test_forked_children_do_not_hold_lifeline_write_ends(self):
        read_end, write_end = corroborate._open_lifeline()
        try:
            self.assertFalse(os.get_inheritable(write_end))
            child = os.fork()
            if child == 0:
                try:
                    os.fstat(write_end)
                except OSError:
                    os._exit(0)
                os._exit(1)
            _, status = os.waitpid(child, 0)
            self.assertEqual(os.waitstatus_to_exitcode(status), 0)
        finally:
            corroborate._close_lifeline(write_end)
            os.close(read_end)
        self.assertEqual(corroborate._LIFELINES, {})

    def _completing_backend(self):
        backend = object.__new__(SystemdScopeBackend)
        backend.systemd_run = "systemd-run"
        backend.systemctl = "systemctl"
        group = self.artifact / "group"
        group.mkdir(exist_ok=True)
        (group / "cgroup.events").write_text("populated 0\n")
        process = mock.Mock()
        process.poll.return_value = 0
        process.communicate.return_value = ("out", "err")
        process.returncode = 0
        return backend, group, process

    def test_off_main_thread_run_leaves_signal_handlers_untouched(self):
        backend, group, process = self._completing_backend()
        before = {sig: signal.getsignal(sig) for sig in (signal.SIGINT, signal.SIGTERM)}
        seen = {}

        def control_group(unit):
            seen[threading.current_thread() is threading.main_thread()] = {
                sig: signal.getsignal(sig) for sig in (signal.SIGINT, signal.SIGTERM)
            }
            return group

        results = []
        with mock.patch(
            "rmgpy.tools.eedf.corroborate.subprocess.Popen", return_value=process
        ), mock.patch.object(backend, "_control_group", side_effect=control_group):
            worker = threading.Thread(
                target=lambda: results.append(backend.run(["solver"], cwd=self.artifact, timeout_s=1.0))
            )
            worker.start()
            worker.join(10.0)
            results.append(backend.run(["solver"], cwd=self.artifact, timeout_s=1.0))
        self.assertEqual([r["cleanup_passed"] for r in results], [True, True])
        # Off the main thread nothing is installed; the main thread's handlers govern.
        self.assertEqual(seen[False], before)
        # On the main thread both are deferred during the run and restored after it.
        self.assertNotEqual(seen[True][signal.SIGINT], before[signal.SIGINT])
        self.assertNotEqual(seen[True][signal.SIGTERM], before[signal.SIGTERM])
        self.assertEqual({sig: signal.getsignal(sig) for sig in before}, before)

    def _fake_systemctl(self):
        # Stand-in for `systemctl --user kill --kill-whom=all --signal=KILL <unit>`:
        # the "unit" is every process whose environment has FAKE_UNIT=<unit>.
        path = self.artifact / "fake-systemctl"
        path.write_text(
            "#!" + sys.executable + "\n"
            "import os, signal, sys\n"
            "unit = sys.argv[-1]\n"
            "with open(os.environ['FAKE_SYSTEMCTL_LOG'], 'a') as log:\n"
            "    log.write(' '.join(sys.argv[1:]) + '\\n')\n"
            "marker = ('FAKE_UNIT=' + unit).encode()\n"
            "for name in os.listdir('/proc'):\n"
            "    if name.isdigit() and int(name) != os.getpid():\n"
            "        try:\n"
            "            with open('/proc/%s/environ' % name, 'rb') as handle:\n"
            "                if marker in handle.read().split(b'\\0'):\n"
            "                    os.kill(int(name), signal.SIGKILL)\n"
            "        except OSError:\n"
            "            pass\n"
        )
        path.chmod(0o755)
        return path

    @staticmethod
    def _alive(pid):
        try:
            state = Path("/proc/%d/stat" % pid).read_text().rsplit(")", 1)[1].split()[0]
        except OSError:
            return False
        return state != "Z"

    def _owner_watch_zombie_owner(self, use_lifeline, c_fork=False):
        # The owner is SIGKILLed and its parent (this test) does not reap it: a
        # zombie owner must not keep the unit alive.
        systemctl = self._fake_systemctl()
        log = self.artifact / "systemctl.log"
        unit = "eedf-corroborate-owner-test-%d-%s.scope" % (os.getpid(), use_lifeline)
        payload = (
            "import os, pathlib, sys, time\n"
            "d = pathlib.Path(sys.argv[1])\n"
            "if os.fork() == 0:\n"
            "    os.setsid()\n"
            "    if os.fork() == 0:\n"
            "        (d / 'descendant').write_text(str(os.getpid()))\n"
            "        time.sleep(60)\n"
            "    os._exit(0)\n"
            "(d / 'payload').write_text(str(os.getpid()))\n"
            "time.sleep(60)\n"
        )
        owner_script = (
            "import os, subprocess, sys, time\n"
            "read_end, write_end = os.pipe() if sys.argv[6] == 'lifeline' else (-1, -1)\n"
            "subprocess.Popen([sys.executable, '-I', '-S', '-c', sys.argv[1], str(os.getpid()), str(read_end),"
            " sys.argv[2], sys.argv[3], '--', sys.executable, '-c', sys.argv[4], sys.argv[5]],"
            " pass_fds=(read_end,) if read_end >= 0 else ())\n"
            "if read_end >= 0:\n"
            "    os.close(read_end)\n"
            "if sys.argv[7] == 'c-fork':\n"
            "    import ctypes\n"
            "    libc = ctypes.CDLL(None)\n"
            "    helper = libc.fork()\n"  # no at-fork hooks: the helper keeps the write end
            "    if helper == 0:\n"
            "        libc.sleep(30)\n"
            "        libc._exit(0)\n"
            "    open(sys.argv[5] + '/helper', 'w').write(str(helper))\n"
            "time.sleep(60)\n"
        )
        env = dict(os.environ, FAKE_UNIT=unit, FAKE_SYSTEMCTL_LOG=str(log))
        owner = subprocess.Popen(
            [sys.executable, "-c", owner_script, _OWNER_WATCH, unit, str(systemctl), payload, str(self.artifact),
             "lifeline" if use_lifeline else "none", "c-fork" if c_fork else "plain"],
            env=env,
        )
        pids = []
        try:
            deadline = time.monotonic() + 10.0
            while time.monotonic() < deadline and not all(
                (self.artifact / name).exists() for name in ("payload", "descendant")
            ):
                time.sleep(0.01)
            pids = [int((self.artifact / name).read_text()) for name in ("payload", "descendant")]
            self.assertTrue(all(self._alive(pid) for pid in pids))
            if c_fork:
                while time.monotonic() < deadline and not (self.artifact / "helper").exists():
                    time.sleep(0.01)
                helper = int((self.artifact / "helper").read_text())
                self.assertTrue(self._alive(helper))
            owner.send_signal(signal.SIGKILL)
            deadline = time.monotonic() + 5.0
            while time.monotonic() < deadline and any(self._alive(pid) for pid in pids):
                time.sleep(0.01)
            self.assertEqual(Path("/proc/%d/stat" % owner.pid).read_text().rsplit(")", 1)[1].split()[0], "Z")
            self.assertEqual([pid for pid in pids if self._alive(pid)], [])
            self.assertIn("--user kill --kill-whom=all --signal=KILL " + unit, log.read_text())
        finally:
            if owner.poll() is None:
                owner.kill()
            owner.wait()
            if c_fork and (self.artifact / "helper").exists():
                pids.append(int((self.artifact / "helper").read_text()))
            for pid in pids:
                if self._alive(pid):
                    os.kill(pid, signal.SIGKILL)

    def test_owner_watch_survives_a_c_level_fork_holding_the_lifeline(self):
        # A helper forked from C (no at-fork hooks) keeps the write end open;
        # the /proc poll still sees the zombie owner die.
        self._owner_watch_zombie_owner(use_lifeline=True, c_fork=True)

    def _run_python(self, code):
        completed = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, timeout=60.0, check=False,
            cwd=self.artifact, env=repository_env(),
        )
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        return completed.stdout

    def test_fork_hook_closes_only_the_lifeline_it_registered(self):
        # An application at-fork hook registered earlier reuses the write end's
        # descriptor number for its own file; the lifeline hook must leave it.
        out = self._run_python(
            "import os\n"
            "state = {}\n"
            "def app_hook():\n"
            "    if 'fd' in state:\n"
            "        reopened = os.open(state['path'], os.O_RDONLY)\n"
            "        os.dup2(reopened, state['fd'])\n"
            "        os.close(reopened)\n"
            "os.register_at_fork(after_in_child=app_hook)\n"
            "from rmgpy.tools.eedf import corroborate\n"
            "unrelated = 0\n"
            "for trial in range(20):\n"
            "    path = 'unrelated-%d' % trial\n"
            "    open(path, 'w').write('app data %d' % trial)\n"
            "    read_end, write_end = corroborate._open_lifeline()\n"
            "    state.update(fd=write_end, path=path)\n"
            "    child = os.fork()\n"
            "    if child == 0:\n"
            "        try:\n"
            "            ok = os.pread(write_end, 64, 0) == ('app data %d' % trial).encode()\n"
            "        except OSError:\n"
            "            ok = False\n"
            "        os._exit(0 if ok else 1)\n"
            "    _, status = os.waitpid(child, 0)\n"
            "    unrelated += os.waitstatus_to_exitcode(status) != 0\n"
            "    state.clear()\n"
            "    corroborate._close_lifeline(write_end)\n"
            "    os.close(read_end)\n"
            "print('closed_unrelated', unrelated)\n"
        )
        self.assertIn("closed_unrelated 0", out)

    def test_foreign_dispositions_are_left_byte_identical(self):
        # A C-level SIGTERM handler installed with ctypes (getsignal() still
        # reports SIG_DFL) and SA_RESTART on SIGINT via siginterrupt(False):
        # after run() returns or raises, both raw sigactions are unchanged,
        # and the C handler was never replaced during the run.
        out = self._run_python(
            "import ctypes, signal, subprocess\n"
            "from pathlib import Path\n"
            "from unittest import mock\n"
            "from rmgpy.tools.eedf import corroborate as c\n"
            "libc = ctypes.CDLL(None)\n"
            "handler = ctypes.CFUNCTYPE(None, ctypes.c_int)(lambda signum: None)\n"
            "class Sigaction(ctypes.Structure):\n"  # glibc's struct sigaction
            "    _fields_ = [('handler', ctypes.c_void_p), ('mask', ctypes.c_ulong * 16),\n"
            "                ('flags', ctypes.c_int), ('restorer', ctypes.c_void_p)]\n"
            "action = Sigaction()\n"
            "action.handler = ctypes.cast(handler, ctypes.c_void_p).value\n"
            "assert libc.sigaction(signal.SIGTERM, ctypes.byref(action), None) == 0\n"
            "signal.siginterrupt(signal.SIGINT, False)\n"
            "def raw(signum):\n"
            "    return bytes(c._raw_sigaction(signum))\n"
            "before = {s: raw(s) for s in (signal.SIGTERM, signal.SIGINT)}\n"
            "assert c._raw_sigaction(signal.SIGINT).flags & 0x10000000\n"
            "group = Path('group'); group.mkdir(exist_ok=True)\n"
            "(group / 'cgroup.events').write_text('populated 0\\n')\n"
            "backend = object.__new__(c.SystemdScopeBackend)\n"
            "backend.systemd_run = 'systemd-run'; backend.systemctl = 'systemctl'\n"
            "during = []\n"
            "def discovery(unit, fail):\n"
            "    during.append(raw(signal.SIGTERM) == before[signal.SIGTERM])\n"
            "    if fail:\n"
            "        raise OSError('discovery failed')\n"
            "    return group\n"
            "process = mock.Mock(); process.poll.return_value = 0\n"
            "process.communicate.return_value = ('', ''); process.returncode = 0\n"
            "unit_kill = mock.Mock(return_value=subprocess.CompletedProcess([], 0, '', ''))\n"
            "for fail in (False, True):\n"
            "    with mock.patch.object(c.subprocess, 'Popen', return_value=process), \\\n"
            "            mock.patch.object(c.subprocess, 'run', unit_kill), \\\n"
            "            mock.patch.object(backend, '_control_group', side_effect=lambda u: discovery(u, fail)):\n"
            "        try:\n"
            "            backend.run(['solver'], cwd='.', timeout_s=1.0)\n"
            "            outcome = 'return'\n"
            "        except OSError:\n"
            "            outcome = 'raise'\n"
            "    after = {s: raw(s) for s in (signal.SIGTERM, signal.SIGINT)}\n"
            "    print(outcome, 'identical', after == before, 'untouched_during', all(during))\n"
        )
        self.assertIn("return identical True untouched_during True", out)
        self.assertIn("raise identical True untouched_during True", out)

    def test_rdkit_sigint_handler_left_by_embedding_is_foreign_and_untouched(self):
        # The cross-test state behind a combined census + corroborate run:
        # RDKit's EmbedMolecule leaves its own SIGINT handler while
        # getsignal() still reports default_int_handler. run() must leave that
        # disposition byte-identical (it is not deferred), and once Python's
        # handler is re-established a SIGINT during cleanup is delivered as
        # KeyboardInterrupt after the unit is verified, not as SIG_DFL.
        out = self._run_python(
            "import ctypes, os, signal, subprocess\n"
            "from pathlib import Path\n"
            "from unittest import mock\n"
            "from rdkit import Chem\n"
            "from rdkit.Chem import AllChem\n"
            "from rmgpy.tools.eedf import corroborate as c\n"
            "python_handler = bytes(c._raw_sigaction(signal.SIGINT))\n"
            "molecule = Chem.AddHs(Chem.MolFromSmiles('CCO'))\n"
            "AllChem.EmbedMolecule(molecule, randomSeed=1)\n"
            "left = bytes(c._raw_sigaction(signal.SIGINT))\n"
            "print('rdkit_left_foreign_handler', left != python_handler,\n"
            "      'stale_view', signal.getsignal(signal.SIGINT) is signal.default_int_handler)\n"
            "group = Path('group'); group.mkdir(exist_ok=True)\n"
            "events = group / 'cgroup.events'\n"
            "backend = object.__new__(c.SystemdScopeBackend)\n"
            "backend.systemd_run = 'systemd-run'; backend.systemctl = 'systemctl'\n"
            "process = mock.Mock(); process.poll.return_value = 0\n"
            "process.communicate.return_value = ('', ''); process.returncode = 0\n"
            "during = []\n"
            "def discovery(unit):\n"
            "    during.append(bytes(c._raw_sigaction(signal.SIGINT)))\n"
            "    events.write_text('populated 1\\n')\n"
            "    return group\n"
            "def unit_kill(argv, **kwargs):\n"
            "    if send:\n"
            "        os.kill(os.getpid(), signal.SIGINT)\n"
            "    events.write_text('populated 0\\n')\n"
            "    return subprocess.CompletedProcess(argv, 0, '', '')\n"
            "def run():\n"
            "    with mock.patch.object(c.subprocess, 'Popen', return_value=process), \\\n"
            "            mock.patch.object(c.subprocess, 'run', side_effect=unit_kill), \\\n"
            "            mock.patch.object(backend, '_control_group', side_effect=discovery):\n"
            "        try:\n"
            "            backend.run(['solver'], cwd='.', timeout_s=0.05)\n"
            "            return 'return'\n"
            "        except KeyboardInterrupt:\n"
            "            return 'KeyboardInterrupt verified=%s' % ('populated 0' in events.read_text())\n"
            "        except c.CorroborationError as exc:\n"
            "            return 'CorroborationError'\n"
            "send = False\n"
            "outcome = run()\n"
            "print('foreign_untouched_during', during[-1] == left,\n"
            "      'identical_after', bytes(c._raw_sigaction(signal.SIGINT)) == left)\n"
            "signal.signal(signal.SIGINT, signal.default_int_handler)\n"
            "send = True\n"
            "print('python_owned', run())\n"
        )
        self.assertIn("rdkit_left_foreign_handler True stale_view True", out)
        self.assertIn("foreign_untouched_during True identical_after True", out)
        self.assertIn("python_owned KeyboardInterrupt verified=True", out)

    def _run_behind_pythons_back(self, raw_handler):
        # Create, in a fresh interpreter, the state RDKit's handler leaves
        # after it has swallowed one SIGINT: the kernel's SIGINT disposition is
        # raw_handler (0 = SIG_DFL, 1 = SIG_IGN) while Python's table still
        # says default_int_handler. The test fixture cannot hide it.
        code = (
            "import os, signal, subprocess, sys\n"
            "from pathlib import Path\n"
            "from unittest import mock\n"
            "from rmgpy.tools.eedf import corroborate as c\n"
            "action = c._raw_sigaction(signal.SIGINT)\n"
            "action.handler = %d\n"
            "action.flags |= 0x10000000\n"  # SA_RESTART, as RDKit's signal() leaves it
            "assert c._rt_sigaction(signal.SIGINT, action, None) == 0\n"
            "assert signal.getsignal(signal.SIGINT) is signal.default_int_handler\n"
            "before = bytes(c._raw_sigaction(signal.SIGINT))\n"
            "group = Path('group'); group.mkdir(exist_ok=True)\n"
            "events = group / 'cgroup.events'\n"
            "verified = Path('verified')\n"
            "backend = object.__new__(c.SystemdScopeBackend)\n"
            "backend.systemd_run = 'systemd-run'; backend.systemctl = 'systemctl'\n"
            "process = mock.Mock(); process.poll.return_value = 0\n"
            "process.communicate.return_value = ('', ''); process.returncode = 0\n"
            "empty = c.SystemdScopeBackend._empty\n"
            "def checked(group):\n"
            "    result = empty(group)\n"
            "    if result:\n"
            "        verified.write_text('unit verified empty')\n"
            "    return result\n"
            "def discovery(unit):\n"
            "    events.write_text('populated 1\\n')\n"
            "    return group\n"
            "def unit_kill(argv, **kwargs):\n"
            "    if send:\n"
            "        os.kill(os.getpid(), signal.SIGINT)\n"
            "    events.write_text('populated 0\\n')\n"
            "    return subprocess.CompletedProcess(argv, 0, '', '')\n"
            "def run():\n"
            "    with mock.patch.object(c.subprocess, 'Popen', return_value=process), \\\n"
            "            mock.patch.object(c.subprocess, 'run', side_effect=unit_kill), \\\n"
            "            mock.patch.object(backend, '_empty', side_effect=checked), \\\n"
            "            mock.patch.object(backend, '_control_group', side_effect=discovery):\n"
            "        try:\n"
            "            backend.run(['solver'], cwd='.', timeout_s=0.05)\n"
            "            return 'return'\n"
            "        except BaseException as exc:\n"
            "            return type(exc).__name__\n"
            "send = False\n"
            "table = signal.getsignal(signal.SIGINT)\n"
            "print('quiet', run(), 'table', 'SIG_DFL' if signal.getsignal(signal.SIGINT) == signal.SIG_DFL else\n"
            "      'python' if signal.getsignal(signal.SIGINT) is signal.default_int_handler else 'other',\n"
            "      'raw_identical', bytes(c._raw_sigaction(signal.SIGINT)) == before, flush=True)\n"
            "verified.unlink()\n"
            "send = True\n"
            "outcome = run()\n"
            "print('survived', outcome, 'raw_identical', bytes(c._raw_sigaction(signal.SIGINT)) == before, flush=True)\n"
        ) % raw_handler
        completed = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, timeout=60.0, check=False,
            cwd=self.artifact, env=repository_env(),
        )
        return completed, (self.artifact / "verified").exists()

    def test_sigint_reset_to_default_behind_pythons_back_is_deferred_as_default(self):
        # The prior disposition really is SIG_DFL: the run is cleaned up and
        # verified first, then SIGINT takes its default action.
        completed, verified = self._run_behind_pythons_back(0)
        # Python's table is left saying what the kernel does.
        self.assertIn("quiet return table SIG_DFL raw_identical True", completed.stdout, completed.stderr)
        self.assertNotIn("survived", completed.stdout)
        self.assertEqual(completed.returncode, -signal.SIGINT, completed.stdout + completed.stderr)
        self.assertTrue(verified, completed.stdout + completed.stderr)

    _NATIVE_SIGINT_PRELUDE = (
        "import ctypes, errno, os, signal, subprocess, sys, threading, time\n"
        "from pathlib import Path\n"
        "from unittest import mock\n"
        "from rmgpy.tools.eedf import corroborate as c\n"
        "libc = ctypes.CDLL(None, use_errno=True)\n"
        "libc.read.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_size_t]\n"
        "libc.read.restype = ctypes.c_ssize_t\n"
        "class Val(ctypes.Union):\n"
        "    _fields_ = [('integer', ctypes.c_int), ('pointer', ctypes.c_void_p)]\n"
        "class Info(ctypes.Structure):\n"
        "    _fields_ = [('signum', ctypes.c_int), ('errno', ctypes.c_int), ('code', ctypes.c_int),\n"
        "                ('padding', ctypes.c_int), ('pid', ctypes.c_int), ('uid', ctypes.c_uint), ('value', Val)]\n"
        "libc.sigqueue.argtypes = [ctypes.c_int, ctypes.c_int, Val]\n"
        "events = []\n"
        "@ctypes.CFUNCTYPE(None, ctypes.c_int)\n"
        "def simple(signum):\n"
        "    events.append(('simple', signum, threading.get_native_id()))\n"
        "@ctypes.CFUNCTYPE(None, ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p)\n"
        "def with_info(signum, info, context):\n"
        "    x = ctypes.cast(info, ctypes.POINTER(Info)).contents\n"
        "    events.append(('siginfo', signum, x.code, x.value.integer))\n"
        # First run in this process: SIGTERM ignored and SIGINT on a C handler,
        # so no deferred signal is at SIG_DFL and CPython's address cannot be learnt.
        "signal.signal(signal.SIGTERM, signal.SIG_IGN)\n"
        "signal.signal(signal.SIGINT, signal.default_int_handler)\n"
        "assert c._PYTHON_TRAMPOLINE is None\n"
        "def install(handler, flags):\n"
        "    action = c._raw_sigaction(signal.SIGINT)\n"
        "    action.handler = ctypes.cast(handler, ctypes.c_void_p).value\n"
        "    action.flags |= flags\n"
        "    assert c._rt_sigaction(signal.SIGINT, action, None) == 0\n"
        "    return bytes(c._raw_sigaction(signal.SIGINT))\n"
        "group = Path('group'); group.mkdir(exist_ok=True)\n"
        "(group / 'cgroup.events').write_text('populated 0\\n')\n"
        "backend = object.__new__(c.SystemdScopeBackend)\n"
        "backend.systemd_run = 'systemd-run'; backend.systemctl = 'systemctl'\n"
        "process = mock.Mock(); process.poll.return_value = 0\n"
        "process.communicate.return_value = ('', ''); process.returncode = 0\n"
        "src, first = __import__('inspect').getsourcelines(c.SystemdScopeBackend.run)\n"
        "unit_line = next(first + i for i, l in enumerate(src) if 'unit = \"eedf-corroborate-' in l)\n"
        "cleanup_src, cleanup_first = __import__('inspect').getsourcelines(c.SystemdScopeBackend._kill_and_verify)\n"
        "cleanup_lines = range(cleanup_first + 1, cleanup_first + len(cleanup_src))\n"
        "def run(on_line):\n"
        "    codes = (c.SystemdScopeBackend.run.__code__, c.SystemdScopeBackend._kill_and_verify.__code__)\n"
        "    def trace(frame, event, arg):\n"
        "        if frame.f_code in codes and event == 'line':\n"
        "            on_line(frame.f_lineno)\n"
        "        return trace\n"
        "    with mock.patch.object(c.subprocess, 'Popen', return_value=process), \\\n"
        "            mock.patch.object(backend, '_control_group', return_value=group):\n"
        "        sys.settrace(trace)\n"
        "        try:\n"
        "            backend.run(['solver'], cwd='.', timeout_s=1.0)\n"
        "        finally:\n"
        "            sys.settrace(None)\n"
    )

    def test_c_handler_without_bootstrap_is_left_foreign_throughout(self):
        # No bootstrap is available, so the C handler stays foreign: its raw
        # action is identical at every executed line of the run.
        out = self._run_python(
            self._NATIVE_SIGINT_PRELUDE
            + "before = install(simple, 0x10000000)\n"
            "changed = []\n"
            "lines = []\n"
            "def check(line):\n"
            "    lines.append(line)\n"
            "    if bytes(c._raw_sigaction(signal.SIGINT)) != before:\n"
            "        changed.append(line)\n"
            "run(check)\n"
            "print('learnt', c._PYTHON_TRAMPOLINE is not None, 'lines', len(set(lines)) > 20,\n"
            "      'changed', changed, 'identical_after', bytes(c._raw_sigaction(signal.SIGINT)) == before)\n"
        )
        self.assertIn("learnt False lines True changed [] identical_after True", out)

    def test_c_handler_keeps_sa_restart_for_native_reads_during_run(self):
        # The r259 bootstrap_eintr probe: a native read in another thread under
        # a foreign SA_RESTART handler is not interrupted by a thread-directed
        # SIGINT, sent wherever the run first alters the action, or else right
        # after ownership is decided.
        out = self._run_python(
            self._NATIVE_SIGINT_PRELUDE
            + "before = install(simple, 0x10000000)\n"
            "read_end, write_end = os.pipe()\n"
            "started = threading.Event(); done = threading.Event(); data = {}\n"
            "def reader():\n"
            "    data['tid'] = threading.get_native_id(); started.set()\n"
            "    buffer = ctypes.create_string_buffer(1)\n"
            "    data['read'] = libc.read(read_end, buffer, 1); data['errno'] = ctypes.get_errno(); done.set()\n"
            "thread = threading.Thread(target=reader); thread.start(); started.wait(); time.sleep(0.02)\n"
            # Inject at the first line where the action differs (none, if it is
            # never touched), else right after ownership is decided.
            "def inject(line):\n"
            "    if 'sent' in data:\n"
            "        return\n"
            "    untouched = bytes(c._raw_sigaction(signal.SIGINT)) == before\n"
            "    if line == unit_line or not untouched:\n"
            "        data['sent'] = untouched\n"
            "        assert libc.syscall(234, os.getpid(), data['tid'], int(signal.SIGINT)) == 0\n"
            "        time.sleep(0.05)\n"
            "        data['returned_early'] = done.is_set()\n"
            "run(inject)\n"
            "os.write(write_end, b'x'); thread.join(5)\n"
            "print('untouched_at_injection', data['sent'], 'returned_early', data['returned_early'],\n"
            "      'read', data['read'], 'errno', data['errno'],\n"
            "      'handled', [e[0] for e in events], 'identical_after', bytes(c._raw_sigaction(signal.SIGINT)) == before)\n"
        )
        self.assertIn(
            "untouched_at_injection True returned_early False read 1 errno 0 handled ['simple'] identical_after True",
            out,
        )

    def test_c_siginfo_handler_keeps_queued_signal_information(self):
        # The r259 siginfo probe: a queued SIGINT during the run reaches the
        # foreign SA_SIGINFO handler with its SI_QUEUE code and value.
        out = self._run_python(
            self._NATIVE_SIGINT_PRELUDE
            + "before = install(with_info, 0x10000000 | 4)\n"
            "def inject(line):\n"
            "    if line == unit_line and not events:\n"
            "        assert libc.sigqueue(os.getpid(), signal.SIGINT, Val(integer=23456)) == 0\n"
            "run(inject)\n"
            "print('events', events, 'identical_after', bytes(c._raw_sigaction(signal.SIGINT)) == before)\n"
        )
        self.assertIn("events [('siginfo', 2, -1, 23456)] identical_after True", out)

    def test_early_python_handlers_on_both_signals_defer_on_first_run_without_a_write(self):
        # R5: both signals sit on Python handlers in a process whose first run
        # this is, so no SIG_DFL install can teach CPython's address. It is
        # identified from the interpreter image, reading only: both raw
        # actions and Python's table are identical at every executed line up
        # to the first defer install, and SIGTERM is then deferred (with
        # SIGINT, which it subsumes) until the unit is verified.
        out = self._run_python(
            self._NATIVE_SIGINT_PRELUDE
            + "term_calls = []\n"
            "signal.signal(signal.SIGTERM, lambda signum, frame: term_calls.append(signum))\n"
            "def state():\n"
            "    return {s: (bytes(c._raw_sigaction(s)), signal.getsignal(s)) for s in (signal.SIGTERM, signal.SIGINT)}\n"
            "before = state()\n"
            "install_line = next(first + i for i, l in enumerate(src)\n"
            "                    if 'signal.signal(guarded_signal, defer_signal)' in l)\n"
            "views = {s: signal.getsignal(s) for s in (signal.SIGTERM, signal.SIGINT)}\n"
            "executed = []; data = {}\n"
            "def check(line):\n"
            "    executed.append((line, state() == before))\n"
            "    if line in cleanup_lines and 'sent' not in data:\n"
            "        data['sent'] = True\n"
            "        os.kill(os.getpid(), signal.SIGTERM)\n"
            "        os.kill(os.getpid(), signal.SIGINT)\n"
            "        data['at_send'] = list(term_calls)\n"
            "        data['owned'] = all(signal.getsignal(s) is not views[s] for s in views)\n"
            "try:\n"
            "    run(check)\n"
            "    outcome = 'return'\n"
            "except KeyboardInterrupt:\n"
            "    outcome = 'KeyboardInterrupt'\n"
            "first_write = next(i for i, (line, same) in enumerate(executed) if not same)\n"
            "print('untouched_until_install', executed[first_write - 1][0] == install_line,\n"
            "      'both_owned', data['owned'], 'deferred', data['at_send'] == [] and term_calls == [signal.SIGTERM],\n"
            "      'outcome', outcome,\n"
            "      'learnt', c._PYTHON_TRAMPOLINE == c._raw_sigaction(signal.SIGINT).handler,\n"
            "      'identical_after', state() == before)\n"
        )
        self.assertIn(
            "untouched_until_install True both_owned True deferred True outcome return learnt True identical_after True",
            out,
        )

    def _foreign_sigint_beside_python_sigterm(self, setup):
        # ``setup`` leaves a non-CPython handler on SIGINT (as ``before``) while
        # Python's table still holds default_int_handler. SIGTERM carries a
        # Python handler, so CPython's address comes from the image path. The
        # SIGINT action is identical at every executed line and after; SIGTERM
        # is deferred.
        return self._run_python(
            self._NATIVE_SIGINT_PRELUDE
            + setup
            + "assert callable(signal.getsignal(signal.SIGINT))\n"
            "term_calls = []\n"
            "signal.signal(signal.SIGTERM, lambda signum, frame: term_calls.append(signum))\n"
            "term_before = bytes(c._raw_sigaction(signal.SIGTERM))\n"
            "changed = []; data = {}\n"
            "sigint_view = signal.getsignal(signal.SIGINT)\n"
            "def check(line):\n"
            "    if bytes(c._raw_sigaction(signal.SIGINT)) != before or signal.getsignal(signal.SIGINT) is not sigint_view:\n"
            "        changed.append(line)\n"
            "    if line in cleanup_lines and 'sent' not in data:\n"
            "        data['sent'] = True\n"
            "        os.kill(os.getpid(), signal.SIGTERM)\n"
            "        data['at_send'] = list(term_calls)\n"
            "run(check)\n"
            "print('sigint_changed', changed,\n"
            "      'sigint_identical_after', bytes(c._raw_sigaction(signal.SIGINT)) == before,\n"
            "      'sigterm_deferred', data['at_send'] == [] and term_calls == [signal.SIGTERM],\n"
            "      'learnt_sigterm_handler', c._PYTHON_TRAMPOLINE == c._raw_sigaction(signal.SIGTERM).handler,\n"
            "      'sigterm_identical_after', bytes(c._raw_sigaction(signal.SIGTERM)) == term_before)\n"
        )

    _FOREIGN_SIGINT_RESULT = (
        "sigint_changed [] sigint_identical_after True sigterm_deferred True "
        "learnt_sigterm_handler True sigterm_identical_after True"
    )

    def test_faulthandler_on_sigint_is_left_foreign(self):
        # faulthandler_user lies in the interpreter image but is not
        # signal_handler.
        out = self._foreign_sigint_beside_python_sigterm(
            "import faulthandler\n"
            "faulthandler.register(signal.SIGINT, file=open('faulthandler.log', 'w'))\n"
            "before = bytes(c._raw_sigaction(signal.SIGINT))\n"
            "print('restart', bool(c._raw_sigaction(signal.SIGINT).flags & 0x10000000))\n"
        )
        self.assertIn("restart True", out)
        self.assertIn(self._FOREIGN_SIGINT_RESULT, out)

    def test_rdkit_sigint_handler_is_left_foreign_beside_a_python_sigterm(self):
        out = self._foreign_sigint_beside_python_sigterm(
            "from rdkit import Chem\n"
            "from rdkit.Chem import AllChem\n"
            "python_sigint = bytes(c._raw_sigaction(signal.SIGINT))\n"
            "AllChem.EmbedMolecule(Chem.AddHs(Chem.MolFromSmiles('CCO')), randomSeed=1)\n"
            "before = bytes(c._raw_sigaction(signal.SIGINT))\n"
            "print('rdkit_handler', before != python_sigint)\n"
        )
        self.assertIn("rdkit_handler True", out)
        self.assertIn(self._FOREIGN_SIGINT_RESULT, out)

    def test_ctypes_handler_with_cpythons_flags_is_left_foreign(self):
        # Exactly PyOS_setsig's flags; it is still not signal_handler.
        out = self._foreign_sigint_beside_python_sigterm(
            "python_flags = c._raw_sigaction(signal.SIGINT).flags\n"
            "before = install(simple, 0)\n"
            "print('same_flags', c._raw_sigaction(signal.SIGINT).flags == python_flags)\n"
        )
        self.assertIn("same_flags True", out)
        self.assertIn(self._FOREIGN_SIGINT_RESULT, out)

    def test_interruptible_faulthandler_on_sigint_is_left_foreign_with_nothing_cached(self):
        # r260: faulthandler.register(chain=False) then siginterrupt(True)
        # leaves faulthandler_user in the interpreter image with exactly
        # PyOS_setsig's flags, and SIGTERM ignored gives no SIG_DFL install.
        # It is still not signal_handler: SIGINT's raw action and table entry
        # are identical at every executed line, and no address is cached. A
        # later run with genuine Python handlers then defers SIGTERM.
        out = self._run_python(
            self._NATIVE_SIGINT_PRELUDE
            + "import faulthandler\n"
            "faulthandler.register(signal.SIGINT, file=open('faulthandler.log', 'w'), chain=False)\n"
            "signal.siginterrupt(signal.SIGINT, True)\n"
            "before = (bytes(c._raw_sigaction(signal.SIGINT)), signal.getsignal(signal.SIGINT))\n"
            "print('pyos_setsig_flags', c._raw_sigaction(signal.SIGINT).flags & ~0x0c000000 == 0)\n"
            "changed = []\n"
            "def check(line):\n"
            "    if (bytes(c._raw_sigaction(signal.SIGINT)), signal.getsignal(signal.SIGINT)) != before:\n"
            "        changed.append(line)\n"
            "run(check)\n"
            "print('changed', changed, 'identical_after',\n"
            "      (bytes(c._raw_sigaction(signal.SIGINT)), signal.getsignal(signal.SIGINT)) == before,\n"
            "      'cached', c._PYTHON_TRAMPOLINE)\n"
            "signal.signal(signal.SIGINT, signal.default_int_handler)\n"
            "term_calls = []\n"
            "signal.signal(signal.SIGTERM, lambda signum, frame: term_calls.append(signum))\n"
            "data = {}\n"
            "def inject(line):\n"
            "    if line in cleanup_lines and 'sent' not in data:\n"
            "        data['sent'] = True\n"
            "        os.kill(os.getpid(), signal.SIGTERM)\n"
            "        data['at_send'] = list(term_calls)\n"
            "run(inject)\n"
            "print('later_deferred', data['at_send'] == [] and term_calls == [signal.SIGTERM],\n"
            "      'later_learnt', c._PYTHON_TRAMPOLINE == c._raw_sigaction(signal.SIGINT).handler)\n"
        )
        self.assertIn("pyos_setsig_flags True", out)
        self.assertIn("changed [] identical_after True cached None", out)
        self.assertIn("later_deferred True later_learnt True", out)

    def test_stripped_interpreter_fails_closed_without_a_write(self):
        # No signal_handler symbol: early Python handlers on both signals are
        # left alone. Neither raw action nor table entry changes at any line.
        out = self._run_python(
            self._NATIVE_SIGINT_PRELUDE
            + "c._elf_function_offset = lambda path, name: None\n"
            "signal.signal(signal.SIGTERM, lambda signum, frame: None)\n"
            "def state():\n"
            "    return {s: (bytes(c._raw_sigaction(s)), signal.getsignal(s)) for s in (signal.SIGTERM, signal.SIGINT)}\n"
            "before = state()\n"
            "changed = []\n"
            "def check(line):\n"
            "    if state() != before:\n"
            "        changed.append(line)\n"
            "run(check)\n"
            "print('changed', changed, 'identical_after', state() == before,\n"
            "      'cached', c._PYTHON_TRAMPOLINE, 'symbol', c._cpython_signal_handler())\n"
        )
        self.assertIn("changed [] identical_after True cached None symbol None", out)

    def test_elf_reader_finds_the_live_signal_handler(self):
        # The pure-Python .symtab reader against the running interpreter: its
        # signal_handler, offset from dladdr's base, is the kernel's SIGINT
        # handler while Python owns SIGINT; faulthandler_user is elsewhere; a
        # stripped copy of the same file yields nothing.
        objcopy = shutil.which("objcopy")
        out = self._run_python(
            "import ctypes, os, shutil, signal, subprocess, sys\n"
            "from rmgpy.tools.eedf import corroborate as c\n"
            "signal.signal(signal.SIGINT, signal.default_int_handler)\n"
            "getsig = ctypes.cast(ctypes.pythonapi.PyOS_getsig, ctypes.c_void_p).value\n"
            "path = c._mapped_file(getsig)\n"
            "offset = c._elf_function_offset(path, 'signal_handler')\n"
            "if offset is None:\n"
            "    print('stripped interpreter'); sys.exit(0)\n"
            "live = c._raw_sigaction(signal.SIGINT).handler\n"
            "print('same_file', os.path.samefile(path, sys.executable),\n"
            "      'live_matches', c._image_base(live) == c._image_base(getsig) and live - c._image_base(live) == offset,\n"
            "      'cached_matches', c._cpython_signal_handler() == live,\n"
            "      'faulthandler_elsewhere', c._elf_function_offset(path, 'faulthandler_user') not in (None, offset))\n"
            "objcopy = %r\n"
            "if objcopy:\n"
            "    subprocess.run([objcopy, '--strip-all', path, 'stripped'], check=True)\n"
            "    print('stripped_copy', c._elf_function_offset('stripped', 'signal_handler'))\n" % objcopy
        )
        if "stripped interpreter" in out:
            self.skipTest("the interpreter has no .symtab")
        self.assertIn("same_file True live_matches True cached_matches True faulthandler_elsewhere True", out)
        if objcopy:
            self.assertIn("stripped_copy None", out)

    def test_address_learnt_from_a_default_install_serves_later_runs(self):
        # The warm-up path: a first run with SIGTERM at SIG_DFL learns
        # CPython's address from its own install; a later run with Python
        # handlers on both signals then defers both without the image path.
        out = self._run_python(
            self._NATIVE_SIGINT_PRELUDE
            + "signal.signal(signal.SIGTERM, signal.SIG_DFL)\n"
            "image = mock.Mock(side_effect=AssertionError('image path used'))\n"
            "c._cpython_handler_by_image = image\n"
            "run(lambda line: None)\n"
            "learnt = c._PYTHON_TRAMPOLINE\n"
            "term_calls = []\n"
            "signal.signal(signal.SIGTERM, lambda signum, frame: term_calls.append(signum))\n"
            "views = {s: signal.getsignal(s) for s in (signal.SIGTERM, signal.SIGINT)}\n"
            "data = {}\n"
            "def check(line):\n"
            "    if line in cleanup_lines and 'sent' not in data:\n"
            "        data['sent'] = True\n"
            "        os.kill(os.getpid(), signal.SIGTERM)\n"
            "        os.kill(os.getpid(), signal.SIGINT)\n"
            "        data['at_send'] = list(term_calls)\n"
            "        data['owned'] = all(signal.getsignal(s) is not views[s] for s in views)\n"
            "try:\n"
            "    run(check)\n"
            "    outcome = 'return'\n"
            "except KeyboardInterrupt:\n"
            "    outcome = 'KeyboardInterrupt'\n"
            "print('learnt', learnt is not None and learnt == c._raw_sigaction(signal.SIGINT).handler,\n"
            "      'both_owned', data['owned'], 'deferred', data['at_send'] == [] and term_calls == [signal.SIGTERM],\n"
            "      'outcome', outcome,\n"
            "      'image_path_calls', image.call_count)\n"
        )
        self.assertIn("learnt True both_owned True deferred True outcome return image_path_calls 0", out)

    def _stale_sigterm_cache(self, cache):
        # SIGTERM's kernel action is SIG_DFL while Python's cached table says
        # None (a C handler was present when the signal module initialised,
        # then native code reset it) or SIG_IGN. The run must defer it as
        # SIG_DFL: verify the unit empty first, then exit -15.
        prelude = (
            "import ctypes, os, signal, subprocess\n"
            "from pathlib import Path\n"
            "from unittest import mock\n"
            "from rmgpy.tools.eedf import corroborate as c\n"
        )
        if cache == "none":
            prelude += (
                "@ctypes.CFUNCTYPE(None, ctypes.c_int)\n"
                "def native(signum):\n"
                "    pass\n"
                "action = c._raw_sigaction(signal.SIGTERM)\n"
                "action.handler = ctypes.cast(native, ctypes.c_void_p).value\n"
                "assert c._rt_sigaction(signal.SIGTERM, action, None) == 0\n"
                # CPython's own initialiser records the C handler as None.
                "init = ctypes.pythonapi.PyInit__signal; init.restype = ctypes.py_object\n"
                "init()\n"
                "assert signal.getsignal(signal.SIGTERM) is None\n"
            )
        else:
            prelude += "signal.signal(signal.SIGTERM, signal.SIG_IGN)\n"
        code = prelude + (
            "action = c._raw_sigaction(signal.SIGTERM); action.handler = 0\n"
            "assert c._rt_sigaction(signal.SIGTERM, action, None) == 0\n"
            "print('cache', signal.getsignal(signal.SIGTERM), 'kernel', c._raw_sigaction(signal.SIGTERM).handler, flush=True)\n"
            "group = Path('group'); group.mkdir(exist_ok=True)\n"
            "(group / 'cgroup.events').write_text('populated 0\\n')\n"
            "backend = object.__new__(c.SystemdScopeBackend)\n"
            "backend.systemd_run = 'systemd-run'; backend.systemctl = 'systemctl'\n"
            "process = mock.Mock(); process.poll.return_value = 0\n"
            "process.communicate.return_value = ('', ''); process.returncode = 0\n"
            "empty = c.SystemdScopeBackend._empty\n"
            "def checked(group):\n"
            "    result = empty(group)\n"
            "    if result:\n"
            "        Path('verified').write_text('yes')\n"
            "    return result\n"
            "def discovery(unit):\n"
            "    os.kill(os.getpid(), signal.SIGTERM)\n"
            "    return group\n"
            "with mock.patch.object(c.subprocess, 'Popen', return_value=process), \\\n"
            "        mock.patch.object(backend, '_control_group', side_effect=discovery), \\\n"
            "        mock.patch.object(backend, '_empty', side_effect=checked):\n"
            "    backend.run(['solver'], cwd='.', timeout_s=1.0)\n"
            "print('survived', flush=True)\n"
        )
        completed = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, timeout=60.0, check=False,
            cwd=self.artifact, env=repository_env(),
        )
        self.assertIn("kernel None", completed.stdout, completed.stderr)
        self.assertEqual(completed.returncode, -signal.SIGTERM, completed.stdout + completed.stderr)
        self.assertNotIn("survived", completed.stdout)
        self.assertTrue((self.artifact / "verified").exists(), completed.stdout + completed.stderr)

    def test_kernel_default_behind_a_cached_none_is_deferred(self):
        self._stale_sigterm_cache("none")

    def test_kernel_default_behind_a_cached_ignore_is_deferred(self):
        self._stale_sigterm_cache("ignore")

    def test_sigint_ignored_behind_pythons_back_stays_ignored(self):
        completed, verified = self._run_behind_pythons_back(1)
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        self.assertIn("quiet return table python raw_identical True", completed.stdout)
        self.assertIn("survived return raw_identical True", completed.stdout)
        self.assertTrue(verified)

    def test_owner_watch_lifeline_kills_unit_while_owner_is_an_unreaped_zombie(self):
        self._owner_watch_zombie_owner(use_lifeline=True)

    def test_owner_watch_proc_fallback_counts_a_zombie_owner_as_dead(self):
        self._owner_watch_zombie_owner(use_lifeline=False)

    def test_owner_watch_execs_payload_with_its_status_and_output(self):
        systemctl = self._fake_systemctl()
        log = self.artifact / "systemctl.log"
        unit = "eedf-corroborate-exec-test-%d.scope" % os.getpid()
        env = dict(os.environ, FAKE_UNIT=unit, FAKE_SYSTEMCTL_LOG=str(log))
        try:
            completed = subprocess.run(
                [sys.executable, "-I", "-S", "-c", _OWNER_WATCH, str(os.getpid()), "-1", unit, str(systemctl),
                 "--", "/bin/sh", "-c", "echo solver-out; echo solver-err >&2; exit 7"],
                capture_output=True,
                text=True,
                timeout=10.0,
                env=env,
            )
            self.assertEqual(completed.returncode, 7)
            self.assertEqual(completed.stdout, "solver-out\n")
            self.assertEqual(completed.stderr, "solver-err\n")
        finally:
            # The watcher outlives the payload until the owner (this process) dies.
            subprocess.run([str(systemctl), "--user", "kill", "--kill-whom=all", "--signal=KILL", unit],
                           env=env, check=True, timeout=10.0)

    def test_owner_watch_refuses_an_owner_that_is_not_its_parent(self):
        systemctl = self._fake_systemctl()
        log = self.artifact / "systemctl.log"
        marker = self.artifact / "payload-ran"
        unit = "eedf-corroborate-orphan-test-%d.scope" % os.getpid()
        stranger = subprocess.Popen(["/bin/sleep", "30"])
        try:
            completed = subprocess.run(
                [sys.executable, "-I", "-S", "-c", _OWNER_WATCH, str(stranger.pid), "-1", unit, str(systemctl),
                 "--", "/bin/sh", "-c", "touch " + str(marker)],
                capture_output=True,
                text=True,
                timeout=10.0,
                env=dict(os.environ, FAKE_UNIT=unit, FAKE_SYSTEMCTL_LOG=str(log)),
            )
            self.assertEqual(completed.returncode, 127)
            deadline = time.monotonic() + 5.0
            while time.monotonic() < deadline and not log.exists():
                time.sleep(0.01)
            self.assertIn("kill --kill-whom=all --signal=KILL " + unit, log.read_text())
            self.assertFalse(marker.exists())
            self.assertIsNone(stranger.poll())
        finally:
            stranger.kill()
            stranger.wait()

    def test_startup_interruption_kills_unowned_launcher(self):
        backend = object.__new__(SystemdScopeBackend)
        backend.systemd_run = "systemd-run"
        backend.systemctl = "systemctl"
        process = mock.Mock()
        process.poll.return_value = None
        with mock.patch(
            "rmgpy.tools.eedf.corroborate.subprocess.Popen", return_value=process
        ), mock.patch.object(
            backend, "_control_group", side_effect=KeyboardInterrupt
        ), mock.patch.object(
            backend, "_kill_and_verify", return_value=False
        ):
            with self.assertRaisesRegex(
                ExecutionBackendUnavailable, "cleanup is unverified"
            ):
                backend.run(["solver"], cwd=self.artifact, timeout_s=1.0)
        process.kill.assert_called_once()
        process.wait.assert_called_once()

    def test_startup_interruption_kills_unknown_unit_after_launcher_exits(self):
        backend = object.__new__(SystemdScopeBackend)
        backend.systemd_run = "systemd-run"
        backend.systemctl = "systemctl"
        process = mock.Mock()
        process.poll.return_value = 0
        unit_kill = mock.Mock(return_value=mock.Mock(returncode=0))
        with mock.patch(
            "rmgpy.tools.eedf.corroborate.subprocess.Popen", return_value=process
        ), mock.patch.object(
            backend, "_control_group", side_effect=KeyboardInterrupt
        ), mock.patch(
            "rmgpy.tools.eedf.corroborate.subprocess.run", unit_kill
        ):
            with self.assertRaisesRegex(
                ExecutionBackendUnavailable, "cleanup is unverified"
            ):
                backend.run(["solver"], cwd=self.artifact, timeout_s=1.0)
        self.assertIn("kill", unit_kill.call_args.args[0])
        self.assertIn("--kill-whom=all", unit_kill.call_args.args[0])

    def test_startup_discovery_failures_kill_unknown_unit(self):
        failures = (
            (None, ExecutionBackendUnavailable),
            (
                ExecutionBackendUnavailable("control-group discovery failed"),
                ExecutionBackendUnavailable,
            ),
            (subprocess.TimeoutExpired(["systemctl", "show"], 2.0), subprocess.TimeoutExpired),
        )
        for discovery, error_type in failures:
            with self.subTest(discovery=discovery):
                backend = object.__new__(SystemdScopeBackend)
                backend.systemd_run = "systemd-run"
                backend.systemctl = "systemctl"
                process = mock.Mock()
                process.poll.return_value = 0
                unit_kill = mock.Mock(return_value=mock.Mock(returncode=0))
                control_group = (
                    mock.patch.object(backend, "_control_group", return_value=None)
                    if discovery is None
                    else mock.patch.object(backend, "_control_group", side_effect=discovery)
                )
                with mock.patch(
                    "rmgpy.tools.eedf.corroborate.subprocess.Popen", return_value=process
                ), control_group, mock.patch(
                    "rmgpy.tools.eedf.corroborate.subprocess.run", unit_kill
                ):
                    with self.assertRaises(error_type):
                        backend.run(["solver"], cwd=self.artifact, timeout_s=1.0)
                self.assertIn("kill", unit_kill.call_args.args[0])
                self.assertIn("--kill-whom=all", unit_kill.call_args.args[0])

    def test_unavailable_cgroup_kill_falls_back_to_whole_unit(self):
        backend = object.__new__(SystemdScopeBackend)
        backend.systemd_run = "systemd-run"
        backend.systemctl = "systemctl"
        group = self.artifact / "solver.scope"
        group.mkdir()
        events = group / "cgroup.events"
        events.write_text("populated 1\n")
        (group / "cgroup.kill").mkdir()
        process = mock.Mock()
        process.poll.return_value = 0
        process.communicate.return_value = ("stdout", "stderr")
        process.returncode = 0

        def kill_unit(*args, **kwargs):
            events.write_text("populated 0\n")
            return mock.Mock(returncode=0, stderr="")

        unit_kill = mock.Mock(side_effect=kill_unit)
        with mock.patch(
            "rmgpy.tools.eedf.corroborate.subprocess.Popen", return_value=process
        ), mock.patch.object(
            backend, "_control_group", return_value=group
        ), mock.patch(
            "rmgpy.tools.eedf.corroborate.subprocess.run", unit_kill
        ):
            result = backend.run(["solver"], cwd=self.artifact, timeout_s=1.0)
        self.assertTrue(result["cleanup_passed"])
        self.assertIn("kill", unit_kill.call_args.args[0])
        self.assertIn("--kill-whom=all", unit_kill.call_args.args[0])

    def test_cleanup_error_still_reaps_live_launcher(self):
        backend = object.__new__(SystemdScopeBackend)
        backend.systemd_run = "systemd-run"
        backend.systemctl = "systemctl"
        process = mock.Mock()
        process.poll.return_value = None
        process.communicate.return_value = ("stdout", "stderr")
        cleanup_timeout = subprocess.TimeoutExpired(["systemctl", "kill"], 2.0)
        with mock.patch(
            "rmgpy.tools.eedf.corroborate.subprocess.Popen", return_value=process
        ), mock.patch.object(
            backend, "_control_group", return_value=self.artifact
        ), mock.patch.object(
            backend, "_kill_and_verify", side_effect=cleanup_timeout
        ):
            with self.assertRaises(subprocess.TimeoutExpired):
                backend.run(["solver"], cwd=self.artifact, timeout_s=1.0)
        process.kill.assert_called_once()
        process.wait.assert_called_once()

    def _run_teardown_interrupt_case(self, *, kill_path, live_launcher):
        root = self.artifact / (kill_path + ("-live" if live_launcher else "-exited"))
        root.mkdir()
        ready = root / "owned-child.pid"
        marker = root / "survived"
        helper = root / "helper.py"
        helper.write_text(
            "import os, pathlib, time\n"
            "if os.fork() == 0:\n"
            "    os.setsid()\n"
            "    if os.fork() == 0:\n"
            "        fd = os.open('/dev/null', os.O_RDWR)\n"
            "        os.dup2(fd, 1)\n"
            "        os.dup2(fd, 2)\n"
            f"        pathlib.Path({str(ready)!r}).write_text(str(os.getpid()))\n"
            "        time.sleep(0.7)\n"
            f"        pathlib.Path({str(marker)!r}).write_text('survived')\n"
            "        os._exit(0)\n"
            "    os._exit(0)\n"
            + ("time.sleep(10)\n" if live_launcher else "")
        )
        group = root / "group"
        group.mkdir()
        events = group / "cgroup.events"
        events.write_text("populated 1\n")
        (group / "cgroup.kill").write_text("")
        backend = object.__new__(SystemdScopeBackend)
        backend.systemd_run = "systemd-run"
        backend.systemctl = "systemctl"
        real_popen = subprocess.Popen
        real_write_text = Path.write_text
        holder = {}
        signals = []
        unrelated = real_popen(["/bin/sleep", "10"])

        def launch(argv, **kwargs):
            holder["process"] = real_popen([sys.executable, str(helper)], **kwargs)
            return holder["process"]

        def control_group(unit):
            deadline = time.monotonic() + 2.0
            while not ready.exists() and time.monotonic() < deadline:
                time.sleep(0.002)
            self.assertTrue(ready.exists())
            if not live_launcher:
                holder["process"].wait(timeout=2.0)
            return group

        def terminate_owned():
            try:
                os.kill(int(ready.read_text()), signal.SIGKILL)
            except ProcessLookupError:
                pass
            if holder["process"].poll() is None:
                holder["process"].kill()
            real_write_text(events, "populated 0\n")

        def write_interface(path, *args, **kwargs):
            if path == group / "cgroup.kill":
                if kill_path == "unit-fallback":
                    raise PermissionError("cgroup.kill unavailable")
                signals.append(signal.SIGINT)
                os.kill(os.getpid(), signal.SIGINT)
                terminate_owned()
            return real_write_text(path, *args, **kwargs)

        def unit_command(argv, **kwargs):
            self.assertIn("--kill-whom=all", argv)
            self.assertIn("--signal=KILL", argv)
            signals.append(signal.SIGINT)
            os.kill(os.getpid(), signal.SIGINT)
            terminate_owned()
            return subprocess.CompletedProcess(argv, 0, "", "")

        try:
            with mock.patch(
                "rmgpy.tools.eedf.corroborate.subprocess.Popen", side_effect=launch
            ), mock.patch.object(
                backend, "_control_group", side_effect=control_group
            ), mock.patch(
                "rmgpy.tools.eedf.corroborate.subprocess.run", side_effect=unit_command
            ), mock.patch.object(Path, "write_text", write_interface):
                with self.assertRaises(KeyboardInterrupt):
                    backend.run(["finite-launcher"], cwd=root, timeout_s=0.05)
            holder["process"].wait(timeout=2.0)
            time.sleep(0.85)
            self.assertEqual(signals, [signal.SIGINT])
            self.assertFalse(marker.exists())
            self.assertIn("populated 0", events.read_text())
            self.assertIsNone(unrelated.poll())
        finally:
            if ready.exists():
                try:
                    os.kill(int(ready.read_text()), signal.SIGKILL)
                except ProcessLookupError:
                    pass
            if holder.get("process") is not None and holder["process"].poll() is None:
                holder["process"].kill()
                holder["process"].wait(timeout=2.0)
            unrelated.terminate()
            unrelated.wait(timeout=2.0)

    def test_cgroup_kill_defers_interrupt_with_exited_launcher(self):
        self._run_teardown_interrupt_case(kill_path="cgroup-kill", live_launcher=False)

    def test_cgroup_kill_defers_interrupt_with_live_launcher(self):
        self._run_teardown_interrupt_case(kill_path="cgroup-kill", live_launcher=True)

    def test_unit_kill_fallback_defers_interrupt_with_exited_launcher(self):
        self._run_teardown_interrupt_case(kill_path="unit-fallback", live_launcher=False)

    def test_unit_kill_fallback_defers_interrupt_with_live_launcher(self):
        self._run_teardown_interrupt_case(kill_path="unit-fallback", live_launcher=True)

    def _run_teardown_class_modes(self, modes):
        root = self.artifact / ("teardown-class-" + modes[0])
        helper = Path(__file__).with_name("teardown_class_helper.py")
        completed = subprocess.run(
            [sys.executable, str(helper), str(root), ",".join(modes)],
            cwd=self.artifact,
            env=repository_env(),
            capture_output=True,
            text=True,
            timeout=10.0 * len(modes),
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr + completed.stdout)
        self.assertEqual(completed.stderr, "")
        for mode in modes:
            with self.subTest(mode=mode):
                case_root = root / mode
                counts = json.loads((case_root / "counts.json").read_text())
                outcome_path = case_root / "outcome.txt"
                outcome = outcome_path.read_text() if outcome_path.exists() else ""
                self.assertFalse((case_root / "survived").exists(), completed.stdout)
                self.assertTrue(counts["verified"], completed.stdout)
                self.assertGreater(counts["cgroup_calls"] + counts["unit_calls"], 0)
                self.assertEqual((case_root / "reviewer.stderr.log").read_text(), "")
                self.assertIn("unrelated_alive True", completed.stdout)
                if mode.startswith("sigterm-"):
                    self.assertIn(f"{mode} exit {-signal.SIGTERM}", completed.stdout)
                elif mode.startswith("sigint-") or mode.startswith("repeated-sigint-"):
                    self.assertTrue(outcome.startswith("KeyboardInterrupt:"), outcome)
        return root, completed

    def test_owned_lifetime_defers_sigterm_until_verified(self):
        self._run_teardown_class_modes(("sigterm-cgroup", "sigterm-unit", "sigterm-verify"))

    def test_owned_lifetime_has_no_post_spawn_interrupt_gap(self):
        # A signal while Popen is still returning (the child already exists) or
        # right after it is recorded and delivered only once the unit is verified.
        self._run_teardown_class_modes(("sigint-inside-spawn", "sigterm-inside-spawn", "sigint-after-spawn"))

    def test_unmasked_idle_thread_cannot_raise_inside_cleanup(self):
        # The signal is directed at an idle thread; Python still runs the
        # deferring handler on the main thread, so cleanup is never interrupted.
        modes = ("idle-thread-sigint-cgroup", "idle-thread-sigterm-cgroup")
        root, completed = self._run_teardown_class_modes(modes)
        self.assertTrue(
            (root / "idle-thread-sigint-cgroup" / "outcome.txt").read_text().startswith("KeyboardInterrupt"),
            completed.stdout,
        )
        self.assertIn(f"idle-thread-sigterm-cgroup exit {-signal.SIGTERM}", completed.stdout)

    def test_teardown_retries_ordinary_kill_and_verify_errors(self):
        modes = ("unit-timeout", "unit-oserror", "verify-oserror")
        self._run_teardown_class_modes(modes)
        expected = {
            "unit-timeout": "TimeoutExpired:",
            "unit-oserror": "OSError:",
            "verify-oserror": "OSError:",
        }
        root = self.artifact / "teardown-class-unit-timeout"
        for mode, prefix in expected.items():
            with self.subTest(mode=mode):
                self.assertTrue((root / mode / "outcome.txt").read_text().startswith(prefix))

    def test_teardown_preserves_prior_interrupt_reproductions(self):
        self._run_teardown_class_modes(
            (
                "sigint-cgroup",
                "repeated-sigint-cgroup",
                "sigint-unit",
                "repeated-sigint-unit",
                "baseexception-cgroup",
                "baseexception-unit",
                "baseexception-verify",
                "sigint-live-launcher",
            )
        )

    def test_stat_errors_retry_unit_kill_and_verify_before_raising(self):
        modes = ("stat-initial-eio", "stat-attempt-eio", "stat-recheck-eio")
        self._run_teardown_class_modes(modes)
        root = self.artifact / "teardown-class-stat-initial-eio"
        for mode in modes:
            with self.subTest(mode=mode):
                case_root = root / mode
                counts = json.loads((case_root / "counts.json").read_text())
                self.assertGreater(counts["unit_calls"], 0)
                self.assertTrue((case_root / "outcome.txt").read_text().startswith("OSError:"))

    def test_real_stat_permission_error_kills_unit_before_raising(self):
        root = self.artifact / "teardown-class-stat-permission"
        helper = Path(__file__).with_name("teardown_class_helper.py")
        completed = subprocess.run(
            [sys.executable, str(helper), str(root), "stat-permission"],
            cwd=self.artifact,
            env=repository_env(),
            capture_output=True,
            text=True,
            timeout=12.0,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr + completed.stdout)
        self.assertEqual(completed.stderr, "")
        case_root = root / "stat-permission"
        counts = json.loads((case_root / "counts.json").read_text())
        self.assertFalse((case_root / "survived").exists(), completed.stdout)
        self.assertGreater(counts["unit_calls"], 0)
        self.assertGreater(counts["verify_calls"], 0)
        self.assertTrue((case_root / "outcome.txt").read_text().startswith("PermissionError:"))
        self.assertIn("unrelated_alive True", completed.stdout)

    def test_mixed_deferred_signals_preserve_termination_and_replay_all(self):
        modes = ("mixed-replay", "mixed-custom-replay")
        root, completed = self._run_teardown_class_modes(modes)
        self.assertIn(f"mixed-replay exit {-signal.SIGTERM}", completed.stdout)
        custom = root / "mixed-custom-replay"
        counts = json.loads((custom / "counts.json").read_text())
        # SIGTERM outranks SIGINT: one replay, SIGINT subsumed.
        self.assertEqual(counts["replayed"], ["SIGTERM"])
        self.assertTrue((custom / "outcome.txt").read_text().startswith("SystemExit: custom SIGTERM replay"))

    def test_sigint_at_release_boundary_does_not_lose_deferred_sigterm(self):
        modes = ("term-int-at-release-term", "term-int-at-release-deliver")
        root, completed = self._run_teardown_class_modes(modes)
        for mode in modes:
            with self.subTest(mode=mode):
                self.assertIn(f"{mode} exit {-signal.SIGTERM}", completed.stdout)
                counts = json.loads((root / mode / "counts.json").read_text())
                self.assertEqual(counts["signals"], ["SIGTERM", "SIGINT"])
                self.assertFalse((root / mode / "survived").exists())
                self.assertFalse((root / mode / "outcome.txt").exists())

    def test_sigint_only_at_release_boundary_is_delivered(self):
        modes = ("int-only-at-release-term", "int-only-at-release-int")
        root, _ = self._run_teardown_class_modes(modes)
        for mode in modes:
            with self.subTest(mode=mode):
                self.assertTrue((root / mode / "outcome.txt").read_text().startswith("KeyboardInterrupt"))
                self.assertFalse((root / mode / "survived").exists())

    def test_subsumption_log_failure_cannot_skip_release(self):
        # A configured FileHandler(delay=True) whose directory is gone or
        # unreadable raises out of logger.info(); the release must still
        # restore both handlers and deliver SIGTERM, leaving nothing pending
        # or blocked.
        modes = (
            "log-failure-control",
            "log-failure-removed-directory",
            "log-failure-permission-denied",
        )
        root, _ = self._run_teardown_class_modes(modes)
        for mode in modes:
            with self.subTest(mode=mode):
                case_root = root / mode
                counts = json.loads((case_root / "counts.json").read_text())
                self.assertEqual(counts["signals"], ["SIGINT", "SIGTERM"])
                self.assertEqual(counts["replayed"], ["SIGTERM"])
                self.assertEqual(counts["pending_after"], [])
                self.assertEqual(counts["blocked_after"], [])
                self.assertTrue(counts["sigint_handler_restored"])
                self.assertTrue(
                    (case_root / "outcome.txt").read_text().startswith("SystemExit: custom SIGTERM replay")
                )


@unittest.skipUnless(
    os.environ.get("EEDF_RUN_CGROUP_TESTS") == "1",
    "requires the deployment host's delegated user systemd scope",
)
class SystemdScopeLifecycleTest(unittest.TestCase):
    """Exercise descendant ownership and cgroup.kill on the deployment host."""

    def setUp(self):
        own_shutdown_signals(self)
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
        result = self.backend.run(["/bin/sh", "-c", "sleep 0.1; exit 0"], cwd=self.root, timeout_s=2)
        self.assertEqual(result["returncode"], 0)
        self.assertTrue(result["cleanup_passed"])

    def test_solver_failure_leaves_no_job_process(self):
        result = self.backend.run(["/bin/sh", "-c", "sleep 0.1; exit 7"], cwd=self.root, timeout_s=2)
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
