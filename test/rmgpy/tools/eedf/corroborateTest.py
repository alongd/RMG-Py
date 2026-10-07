import hashlib
import json
import math
from pathlib import Path
import stat
import tempfile
import time
import unittest
from unittest import mock

import h5py
import numpy as np
from scipy.constants import electron_mass

from rmgpy.tools.eedf.corroborate import (
    CorroborationError,
    _floors,
    _loki_quadrature,
    _quadrature,
    corroborate_row,
)


class CorroborateTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.artifact = Path(self.tmp.name)
        mass = electron_mass / 1e-5
        self._write_collision(
            "elastic.txt",
            "e + Ar(1S0) -> e + Ar(1S0), Elastic",
            "Elastic",
            0.0,
            ((0.0, 1.0), (1.0, 1.0), (2.0, 1.0)),
        )
        self._write_collision(
            "excitation.txt",
            "e + Ar(1S0) -> e + Ar(2P1), Excitation",
            "Excitation",
            0.5,
            ((0.0, 0.0), (1.0, 0.5), (2.0, 0.5)),
        )
        mass_file = self.artifact / "masses.txt"
        mass_file.write_text(f"Ar {mass:.17g}\n")
        self.channels = [
            self._channel(
                "e + Ar(1S0) -> e + Ar(1S0), Elastic",
                "elastic",
                0.0,
                [0.0, 1.0, 2.0],
                [1.0, 1.0, 1.0],
                mass_ratio=1e-5,
            ),
            self._channel(
                "e + Ar(1S0) -> e + Ar(2P1), Excitation",
                "excitation",
                0.5,
                [0.0, 1.0, 2.0],
                [0.0, 0.5, 0.5],
            ),
        ]
        self.processes = [
            dict(
                channel,
                energy=channel["cross_section"]["energy_eV"],
                sigma=channel["cross_section"]["sigma_m2"],
                threshold=channel["threshold_eV"],
            )
            for channel in self.channels
        ]
        input_files = {}
        for name, kind in (
            ("elastic.txt", "cross_section"),
            ("excitation.txt", "cross_section"),
            ("masses.txt", "property"),
        ):
            path = self.artifact / name
            input_files[name] = {
                "path": str(path),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "kind": kind,
            }
        self.spec = {
            "arm": {"gases": ["Ar"], "additive": None, "feed_partial_pressure_Pa": 0.0},
            "Tg_K": 300.0,
            "P_Pa": 100.0,
            "input_files": input_files,
            "envelopes": {},
            "gas_properties": {"mass": "masses.txt", "fraction": ["Ar = 1.0"]},
            "state_properties": {"population": ["Ar(1S0) = 1.0"]},
            "working_conditions": {"excitationFrequency": 0},
            "solver_options": {
                "eedfType": "boltzmann",
                "ionizationOperatorType": "usingSDCS",
                "growthModelType": "temporal",
                "includeEECollisions": False,
                "LXCatFiles": ["elastic.txt", "excitation.txt"],
                "numerics": {"energyGrid": {"maxEnergy": 2.0, "cellNumber": 2}},
            },
            "timeout_s": 1.0,
        }
        energy = np.array([0.5, 1.5])
        edges = np.array([0.0, 1.0, 2.0])
        f0 = np.exp(-energy)
        f0 /= np.sum(np.sqrt(energy) * f0 * np.diff(edges))
        self.row_for_quadrature = {
            "energy_eV": energy,
            "energy_edges_eV": edges,
            "gas_temperature_K": 300.0,
            "target_fractions": np.array([1.0, 1.0]),
            "product_fractions": np.array([0.0, 0.0]),
        }
        baseline = _quadrature(energy, f0, self.processes, self.row_for_quadrature)
        self.manifest = {
            "axes": {"u": [math.log(10.0)]},
            "energy_eV": energy.tolist(),
            "energy_edges_eV": edges.tolist(),
            "branches": ["branch_0"],
            "channel_map": self.channels,
            "floors": {"absolute_power_share": 1e-4},
            "row_inputs": {"Tg_K": 300.0, "P_Pa": 100.0},
        }
        self._write_metadata()
        with h5py.File(self.artifact / "table.h5", "w") as h5:
            group = h5.create_group("branches").create_group("branch_0")
            group.create_dataset("f0", data=[f0])
            group.create_dataset("converged", data=[True])
            group.create_dataset("target_fractions", data=[[1.0, 1.0]])
            group.create_dataset("product_fractions", data=[[0.0, 0.0]])
            group.create_dataset("rate_floors", data=[[1e-40, 1e-40]])
            group.create_dataset("below_floor", data=[[[False, False], [False, False]]])
            group.create_dataset("k_ine", data=[baseline["rates"]])
            group.create_dataset("channel_power", data=[baseline["channel_power"]])
            swarm = group.create_group("swarm")
            swarm.create_dataset("mean_energy_eV", data=[baseline["mean_energy_eV"]])
            swarm.create_dataset("mobility_N", data=[baseline["mobility_N"]])
            power = group.create_group("power_groups")
            power.create_dataset("field", data=[max(abs(baseline["total_power"]), 1.0)])
        self.fake_count = 0

    def tearDown(self):
        self.tmp.cleanup()

    def _write_collision(self, name, description, kind, threshold, values):
        rows = "\n".join(f"{energy} {sigma}" for energy, sigma in values)
        (self.artifact / name).write_text(
            f"PARAM.: E = {threshold} eV\n[{description}]\n----------------\n{rows}\n----------------\n"
        )

    @staticmethod
    def _channel(description, kind, threshold, energy, sigma, mass_ratio=None):
        channel = {
            "description": description,
            "kind": kind,
            "classification": "B",
            "threshold_eV": threshold,
            "target_fraction": 1.0,
            "product_fraction": 0.0,
            "sigma_max_m2": max(sigma),
            "flux_group": kind,
            "reaction": None,
            "cross_section": {"energy_eV": energy, "sigma_m2": sigma},
        }
        if mass_ratio is not None:
            channel["mass_ratio"] = mass_ratio
        return channel

    def _write_metadata(self):
        (self.artifact / "generation_spec.json").write_text(json.dumps(self.spec))
        (self.artifact / "manifest.json").write_text(json.dumps(self.manifest))

    def _add_ionization(self):
        name = "ionization.txt"
        description = "e + Ar(1S0) -> e + e + Ar(+,gnd), Ionization"
        self._write_collision(
            name,
            description,
            "Ionization",
            0.5,
            ((0.0, 0.0), (1.0, 0.25), (2.0, 0.25)),
        )
        path = self.artifact / name
        self.spec["input_files"][name] = {
            "path": str(path),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "kind": "cross_section",
        }
        self.spec["solver_options"]["LXCatFiles"].append(name)
        channel = self._channel(
            description,
            "ionization",
            0.5,
            [0.0, 1.0, 2.0],
            [0.0, 0.25, 0.25],
        )
        channel["opb_eV"] = 0.5
        self.channels.append(channel)
        self.processes.append(
            dict(
                channel,
                energy=channel["cross_section"]["energy_eV"],
                sigma=channel["cross_section"]["sigma_m2"],
                threshold=channel["threshold_eV"],
            )
        )
        self.manifest["channel_map"] = self.channels
        self._write_metadata()
        with h5py.File(self.artifact / "table.h5", "r+") as h5:
            group = h5["branches/branch_0"]
            f0 = np.asarray(group["f0"][0])
            baseline = _quadrature(
                self.row_for_quadrature["energy_eV"],
                f0,
                self.processes,
                dict(
                    self.row_for_quadrature,
                    target_fractions=np.ones(len(self.processes)),
                    product_fractions=np.zeros(len(self.processes)),
                ),
            )
            replacements = {
                "target_fractions": [[1.0] * len(self.processes)],
                "product_fractions": [[0.0] * len(self.processes)],
                "rate_floors": [[1e-40] * len(self.processes)],
                "below_floor": [[[False, False] for _ in range(len(self.processes))]],
                "k_ine": [baseline["rates"]],
                "channel_power": [baseline["channel_power"]],
            }
            for key, value in replacements.items():
                del group[key]
                group.create_dataset(key, data=value)
            group["swarm/mean_energy_eV"][0] = baseline["mean_energy_eV"]
            group["swarm/mobility_N"][0] = baseline["mobility_N"]

    @staticmethod
    def _write_executable(path, body):
        path.write_text("#!/usr/bin/env python3\n" + body)
        path.chmod(path.stat().st_mode | stat.S_IXUSR)

    def _f0(self):
        energy = np.array([0.5, 1.5])
        f0 = np.exp(-energy)
        return f0 / np.sum(np.sqrt(energy) * f0)

    def solver(
        self,
        f0=None,
        name="fake",
        energy=None,
        status="converged",
        expected_config=None,
        **overrides,
    ):
        self.fake_count += 1
        interpreter = self.artifact / f"fake-{name}-{self.fake_count}.py"
        result = {
            "status": status,
            "energy": energy or [0.5, 1.5],
            "f0": (f0 if f0 is not None else self._f0()).tolist(),
            "version": "fake",
        }
        body = "import json\nimport sys\n\nrequest = json.load(open(sys.argv[-1]))\n"
        for key, value in (expected_config or {}).items():
            body += f"assert request[{key!r}] == {value!r}\n"
        body += f"print({json.dumps(result)!r})\n"
        self._write_executable(interpreter, body)
        spec = {"interpreter": str(interpreter), "kind": "bolos"}
        spec.update(overrides)
        return {name: spec}

    def _run_fake_references(self, solvers, node=None):
        """Exercise comparison policy without adding a production solver bypass."""
        with mock.patch(
            "rmgpy.tools.eedf.corroborate._unavailable_reason", return_value=None
        ):
            return corroborate_row(
                self.artifact,
                "branch_0",
                node or {},
                math.log(10.0),
                solvers=solvers,
            )

    def test_agree_and_quadrature_consistency_are_reported_separately(self):
        result = self._run_fake_references(self.solver())
        self.assertEqual(result["overall_verdict"], "AGREE")
        self.assertEqual(result["loki_consistency"]["verdict"], "AGREE")
        self.assertAlmostEqual(
            result["solvers"]["fake"]["ratios"]["mean_energy_eV"], 1.0
        )

    def test_ionization_free_setup_is_admitted_through_production_path(self):
        result = corroborate_row(
            self.artifact,
            "branch_0",
            {},
            math.log(10.0),
            solvers=self.solver(name="bolos"),
        )
        self.assertEqual(result["solvers"]["bolos"]["verdict"], "AGREE")
        self.assertEqual(result["loki_consistency"]["verdict"], "AGREE")
        self.assertEqual(result["overall_verdict"], "AGREE")

    def test_inert_explicit_state_metadata_does_not_refuse_bolos(self):
        self.spec["state_properties"] = {
            "population": ["Ar(1S0) = 1.0", "Ar(2P1) = 0.0"],
            "energy": ["Ar(1S0) = 0.0"],
            "statisticalWeight": ["Ar(1S0) = 1.0", "Ar(2P1) = 3.0"],
        }
        self._write_metadata()
        result = corroborate_row(
            self.artifact,
            "branch_0",
            {},
            math.log(10.0),
            solvers=self.solver(name="bolos"),
        )
        self.assertEqual(result["solvers"]["bolos"]["verdict"], "AGREE")

    def test_equal_sharing_ionization_is_admitted_for_bolos(self):
        self._add_ionization()
        self.spec["solver_options"]["ionizationOperatorType"] = "equalSharing"
        self._write_metadata()
        result = corroborate_row(
            self.artifact,
            "branch_0",
            {},
            math.log(10.0),
            solvers=self.solver(name="bolos"),
        )
        self.assertEqual(result["solvers"]["bolos"]["verdict"], "AGREE")
        self.assertEqual(result["loki_consistency"]["verdict"], "AGREE")

    def test_disagree(self):
        result = self._run_fake_references(self.solver(np.array([10.0, 1.0])))
        self.assertEqual(result["overall_verdict"], "DISAGREE")
        self.assertEqual(
            result["solvers"]["fake"]["quantities"]["rate[0]"]["status"], "DISAGREE"
        )

    def test_materially_negative_reference_eedf_is_unavailable(self):
        result = corroborate_row(
            self.artifact,
            "branch_0",
            {},
            math.log(10.0),
            solvers=self.solver(f0=np.array([1.0, -100.0]), name="bolos"),
        )
        self.assertEqual(result["solvers"]["bolos"]["verdict"], "UNAVAILABLE")
        self.assertIn("negative mass", result["solvers"]["bolos"]["reason"])

    def test_quadrature_refuses_missing_stored_inputs(self):
        complete = dict(
            self.row_for_quadrature,
            target_fractions=np.array([1.0, 1.0]),
            product_fractions=np.array([0.0, 0.0]),
        )
        for key in ("target_fractions", "product_fractions", "gas_temperature_K"):
            with self.subTest(key=key):
                row = dict(complete)
                del row[key]
                with self.assertRaisesRegex(CorroborationError, key):
                    _quadrature(row["energy_eV"], self._f0(), self.processes, row)

    def test_collision_process_comparison_fields_are_required(self):
        row = dict(
            self.row_for_quadrature,
            target_fractions=np.array([1.0]),
            product_fractions=np.array([0.0]),
        )
        for key in ("energy", "sigma", "threshold"):
            with self.subTest(key=key):
                process = dict(self.processes[0])
                del process[key]
                with self.assertRaisesRegex(CorroborationError, key):
                    _quadrature(row["energy_eV"], self._f0(), [process], row)

    def test_all_floor_metadata_must_be_explicit(self):
        complete_row = {
            "rate_floors": np.array([1e-40]),
            "power_groups": {"field": 1.0},
        }
        complete_manifest = {"floors": {"absolute_power_share": 1e-4}}
        cases = (
            ("rate_floors", {}, complete_manifest),
            ("power_groups", {"rate_floors": np.array([1e-40])}, complete_manifest),
            ("floors", complete_row, {}),
        )
        for key, row, manifest in cases:
            with self.subTest(key=key):
                with self.assertRaisesRegex(CorroborationError, key):
                    _floors(row, manifest, 1)

    def test_missing_floor_metadata_refuses(self):
        with h5py.File(self.artifact / "table.h5", "r+") as h5:
            del h5["branches/branch_0/rate_floors"]
        with self.assertRaisesRegex(CorroborationError, "rate_floors"):
            corroborate_row(
                self.artifact,
                "branch_0",
                {},
                math.log(10.0),
                solvers=self.solver(name="bolos"),
            )

    def test_missing_input_file_registry_refuses_by_name(self):
        del self.spec["input_files"]
        self._write_metadata()
        with self.assertRaisesRegex(CorroborationError, "input_files"):
            corroborate_row(
                self.artifact,
                "branch_0",
                {},
                math.log(10.0),
                solvers=self.solver(name="bolos"),
            )

    def test_missing_pinned_temperature_refuses_on_production_path(self):
        del self.spec["Tg_K"]
        del self.manifest["row_inputs"]["Tg_K"]
        self._write_metadata()
        with self.assertRaisesRegex(CorroborationError, "Tg_K"):
            corroborate_row(
                self.artifact,
                "branch_0",
                {},
                math.log(10.0),
                solvers=self.solver(name="bolos"),
            )

    def test_no_averaging(self):
        solvers = self.solver()
        solvers.update(self.solver(np.array([10.0, 1.0]), name="bad"))
        result = self._run_fake_references(solvers)
        self.assertEqual(result["solvers"]["fake"]["verdict"], "AGREE")
        self.assertEqual(result["solvers"]["bad"]["verdict"], "DISAGREE")
        self.assertEqual(result["overall_verdict"], "DISAGREE")

    def test_manifest_collision_arrays_cannot_replace_pinned_bytes(self):
        self.manifest["channel_map"][1]["cross_section"]["sigma_m2"][1] *= 2.0
        self.manifest["channel_map"][1]["sigma_max_m2"] *= 2.0
        self._write_metadata()
        with self.assertRaisesRegex(CorroborationError, "hashed LoKI input"):
            corroborate_row(
                self.artifact, "branch_0", {}, math.log(10.0), solvers=self.solver()
            )

    def test_all_pinned_collision_files_are_loaded(self):
        self.spec["input_files"].pop("excitation.txt")
        self.spec["solver_options"]["LXCatFiles"].remove("excitation.txt")
        self._write_metadata()
        with self.assertRaisesRegex(
            CorroborationError, "physical collision identities"
        ):
            corroborate_row(
                self.artifact, "branch_0", {}, math.log(10.0), solvers=self.solver()
            )

    def test_cantera_is_unavailable_when_grid_is_not_controllable(self):
        result = corroborate_row(
            self.artifact,
            "branch_0",
            {},
            math.log(10.0),
            solvers=self.solver(name="cantera", kind="cantera"),
        )
        self.assertEqual(result["solvers"]["cantera"]["verdict"], "UNAVAILABLE")
        self.assertIn("grid not controllable", result["solvers"]["cantera"]["reason"])

    def test_returned_grid_must_match_the_stored_grid(self):
        solver = self.solver(energy=[0.4, 1.4])
        result = self._run_fake_references(solver)
        self.assertEqual(result["solvers"]["fake"]["verdict"], "UNAVAILABLE")
        self.assertIn("grid", result["solvers"]["fake"]["reason"])

    def test_mixture_is_not_silently_replaced_with_pure_argon(self):
        self.spec["arm"] = {
            "gases": ["Ar", "N2"],
            "additive": "N2",
            "feed_partial_pressure_Pa": 1.0,
        }
        self._write_metadata()
        result = corroborate_row(
            self.artifact,
            "branch_0",
            {},
            math.log(10.0),
            solvers=self.solver(name="bolos", kind="bolos"),
        )
        self.assertIn("pure Ar", result["solvers"]["bolos"]["reason"])

    def test_population_is_not_silently_replaced(self):
        self.spec["state_properties"]["population"] = ["Ar(1S0) = 0.5", "Ar(2P1) = 0.5"]
        self._write_metadata()
        result = corroborate_row(
            self.artifact,
            "branch_0",
            {},
            math.log(10.0),
            solvers=self.solver(name="bolos", kind="bolos"),
        )
        self.assertIn("population", result["solvers"]["bolos"]["reason"])

    def test_nonzero_state_energies_are_not_silently_ignored(self):
        self.spec["state_properties"]["energy"] = ["Ar(1S0) = 0.1"]
        self._write_metadata()
        result = corroborate_row(
            self.artifact,
            "branch_0",
            {},
            math.log(10.0),
            solvers=self.solver(name="bolos", kind="bolos"),
        )
        self.assertIn("state energies", result["solvers"]["bolos"]["reason"])

    def test_electron_electron_collisions_are_not_silently_ignored(self):
        self.spec["solver_options"]["includeEECollisions"] = True
        self._write_metadata()
        result = corroborate_row(
            self.artifact,
            "branch_0",
            {},
            math.log(10.0),
            solvers=self.solver(name="bolos", kind="bolos"),
        )
        self.assertIn("electron-electron", result["solvers"]["bolos"]["reason"])

    def test_loki_operator_is_not_silently_substituted(self):
        self._add_ionization()
        self._write_metadata()
        result = corroborate_row(
            self.artifact,
            "branch_0",
            {},
            math.log(10.0),
            solvers=self.solver(name="bolos", kind="bolos"),
        )
        self.assertIn("usingSDCS", result["solvers"]["bolos"]["reason"])

    def test_unknown_ionization_operator_is_not_accepted(self):
        self._add_ionization()
        self.spec["solver_options"]["ionizationOperatorType"] = "unknown"
        self._write_metadata()
        result = corroborate_row(
            self.artifact,
            "branch_0",
            {},
            math.log(10.0),
            solvers=self.solver(name="bolos", kind="bolos"),
        )
        self.assertIn(
            "not a validated exact match", result["solvers"]["bolos"]["reason"]
        )

    def test_node_temperature_and_pressure_are_passed_to_reference(self):
        self.manifest["axes"] = {
            "u": [math.log(10.0)],
            "Tg_K": [600.0],
            "P_Pa": [200.0],
        }
        self._write_metadata()
        with h5py.File(self.artifact / "table.h5", "r+") as h5:
            datasets = []
            h5.visititems(
                lambda name, value: (
                    datasets.append(name) if isinstance(value, h5py.Dataset) else None
                )
            )
            for path in datasets:
                data = np.asarray(h5[path])
                del h5[path]
                parent, name = path.rsplit("/", 1)
                h5[parent].create_dataset(
                    name, data=data.reshape((1, 1, 1) + data.shape[1:])
                )
        result = self._run_fake_references(
            self.solver(expected_config={"Tg": 600.0, "P": 200.0}),
            node={"Tg_K": 600.0, "P_Pa": 200.0},
        )
        self.assertNotEqual(result["solvers"]["fake"]["verdict"], "UNAVAILABLE")

    def test_envelope_reference_conditions_are_passed_to_reference(self):
        self.spec["envelopes"] = {
            "Tg_K": {"reference": 450.0},
            "P_Pa": {"reference": 150.0},
        }
        self._write_metadata()
        result = self._run_fake_references(
            self.solver(expected_config={"Tg": 450.0, "P": 150.0})
        )
        self.assertNotEqual(result["solvers"]["fake"]["verdict"], "UNAVAILABLE")

    def test_reference_is_compared_to_stored_row_not_recomputed_loki(self):
        with h5py.File(self.artifact / "table.h5", "r+") as h5:
            h5["branches/branch_0/swarm/mean_energy_eV"][0] *= 2.0
        result = self._run_fake_references(self.solver())
        comparison = result["solvers"]["fake"]["quantities"]["mean_energy_eV"]
        self.assertAlmostEqual(comparison["ratio"], 0.5)
        self.assertEqual(comparison["status"], "DISAGREE")
        self.assertEqual(result["loki_consistency"]["verdict"], "DISAGREE")

    def test_unconverged_stored_loki_row_refuses(self):
        with h5py.File(self.artifact / "table.h5", "r+") as h5:
            h5["branches/branch_0/converged"][0] = False
        with self.assertRaisesRegex(
            CorroborationError, "stored LoKI row is not converged"
        ):
            corroborate_row(
                self.artifact, "branch_0", {}, math.log(10.0), solvers=self.solver()
            )

    def test_below_floor_requires_stored_and_every_reference_below(self):
        high_f0 = np.array([1.0, 10.0])
        high = _quadrature(
            np.array([0.5, 1.5]), high_f0, self.processes, self.row_for_quadrature
        )
        with h5py.File(self.artifact / "table.h5", "r+") as h5:
            group = h5["branches/branch_0"]
            stored_rate = abs(group["k_ine"][0, 1])
            group["rate_floors"][0, 1] = (stored_rate + abs(high["rates"][1])) / 2.0
            group["below_floor"][0, 1, 0] = True
        solvers = self.solver(name="low")
        solvers.update(self.solver(high_f0, name="high"))
        result = self._run_fake_references(solvers)
        self.assertNotEqual(
            result["solvers"]["low"]["quantities"]["rate[1]"]["status"], "BELOW_FLOOR"
        )
        self.assertNotEqual(
            result["solvers"]["high"]["quantities"]["rate[1]"]["status"], "BELOW_FLOOR"
        )

    def test_timeout_kills_setsid_descendants_and_keeps_other_results(self):
        marker = self.artifact / "orphan-ran"
        ready = self.artifact / "orphan-ready"
        slow = self.artifact / "slow.py"
        child = (
            "import os\nimport pathlib\nimport time\n\nos.setsid()\n"
            "pathlib.Path(%r).write_text('ready')\ntime.sleep(1.5)\n"
            "pathlib.Path(%r).write_text('orphan')" % (str(ready), str(marker))
        )
        self._write_executable(
            slow,
            "import pathlib\n"
            "import subprocess\nimport sys\nimport time\n\n"
            f"subprocess.Popen([sys.executable, '-c', {child!r}])\n"
            f"while not pathlib.Path({str(ready)!r}).exists():\n"
            "    time.sleep(.005)\n"
            "time.sleep(5)\n",
        )
        solvers = self.solver(name="good")
        solvers["slow"] = dict(
            self.solver(name="slow")["slow"], interpreter=str(slow), timeout_s=0.5
        )
        started = time.monotonic()
        result = self._run_fake_references(solvers)
        elapsed = time.monotonic() - started
        time.sleep(1.6)
        self.assertLess(elapsed, 1.0)
        self.assertTrue(ready.exists())
        self.assertFalse(marker.exists())
        self.assertEqual(result["solvers"]["good"]["verdict"], "AGREE")
        self.assertEqual(result["solvers"]["slow"]["verdict"], "UNAVAILABLE")
        self.assertIn("timeout", result["solvers"]["slow"]["reason"])

    def test_missing_crashed_and_malformed_solvers_are_independently_unavailable(self):
        crash = self.artifact / "crash.py"
        malformed = self.artifact / "malformed.py"
        self._write_executable(crash, "raise SystemExit(3)\n")
        self._write_executable(malformed, "print('not json')\n")
        solvers = self.solver(name="good")
        solvers["missing"] = dict(
            self.solver(name="missing")["missing"],
            interpreter=str(self.artifact / "missing"),
        )
        solvers["crash"] = dict(
            self.solver(name="crash")["crash"], interpreter=str(crash)
        )
        solvers["malformed"] = dict(
            self.solver(name="malformed")["malformed"], interpreter=str(malformed)
        )
        result = self._run_fake_references(solvers)
        self.assertEqual(result["solvers"]["good"]["verdict"], "AGREE")
        for name in ("missing", "crash", "malformed"):
            self.assertEqual(result["solvers"][name]["verdict"], "UNAVAILABLE")

    def test_elastic_power_is_net_transfer_including_thermal_gain(self):
        energy = np.array([0.5, 1.5])
        row = {
            "energy_eV": energy,
            "energy_edges_eV": np.array([0.0, 1.0, 2.0]),
            "target_fractions": np.array([1.0]),
            "product_fractions": np.array([0.0]),
        }
        cold = _quadrature(
            energy, self._f0(), [self.processes[0]], dict(row, gas_temperature_K=0.0)
        )["channel_power"][0]
        hot = _quadrature(
            energy, self._f0(), [self.processes[0]], dict(row, gas_temperature_K=3000.0)
        )["channel_power"][0]
        self.assertLess(hot, cold)
        self.assertNotAlmostEqual(hot, cold, places=8)

    def test_sha_mismatch(self):
        self.spec["input_files"]["elastic.txt"]["sha256"] = "0" * 64
        self._write_metadata()
        with self.assertRaisesRegex(CorroborationError, "SHA-256 mismatch"):
            corroborate_row(
                self.artifact, "branch_0", {}, math.log(10.0), solvers=self.solver()
            )

    def test_superelastic_set_refuses(self):
        description = "e + Ar(1S0) <-> e + Ar(2P1), Excitation"
        self._write_collision(
            "excitation.txt",
            description,
            "Excitation",
            0.5,
            ((0.0, 0.0), (1.0, 0.5), (2.0, 0.5)),
        )
        self.spec["input_files"]["excitation.txt"]["sha256"] = hashlib.sha256(
            (self.artifact / "excitation.txt").read_bytes()
        ).hexdigest()
        self.spec["state_properties"]["statisticalWeight"] = [
            "Ar(1S0) = 1.0",
            "Ar(2P1) = 1.0",
        ]
        self.manifest["channel_map"][1]["description"] = description
        self.manifest["channel_map"][1]["statistical_weight_ratio"] = 1.0
        self._write_metadata()
        result = corroborate_row(
            self.artifact, "branch_0", {}, math.log(10.0), solvers=self.solver()
        )
        self.assertEqual(result["solvers"]["fake"]["verdict"], "UNAVAILABLE")
        self.assertIn("superelastic", result["solvers"]["fake"]["reason"])

    def test_unconverged_reference_is_unavailable(self):
        solver = self.solver(status="unconverged")
        result = self._run_fake_references(solver)
        self.assertEqual(result["overall_verdict"], "UNAVAILABLE")

    def test_solver_process_is_isolated_from_caller(self):
        isolated = self.artifact / "isolated.py"
        self._write_executable(
            isolated,
            "import os\nimport signal\n\nos.killpg(os.getpgrp(), signal.SIGTERM)\n",
        )
        solver = self.solver()
        solver["fake"]["interpreter"] = str(isolated)
        result = self._run_fake_references(solver)
        self.assertEqual(result["solvers"]["fake"]["verdict"], "UNAVAILABLE")

    def test_loki_consistency_uses_sdcs_ionization_power(self):
        process = {
            "description": "e + Ar(1S0) -> e + e + Ar+, Ionization",
            "kind": "ionization",
            "product": "Ar+",
            "threshold": 0.5,
            "energy": [0.0, 1.0, 2.0],
            "sigma": [0.0, 0.5, 0.5],
            "opb_eV": 0.25,
        }
        row = dict(
            self.row_for_quadrature,
            f0=self._f0(),
            target_fractions=np.array([1.0]),
            product_fractions=np.array([0.0]),
        )
        sdcs = _loki_quadrature(row, [process])["channel_power"][0]
        threshold_loss = _quadrature(row["energy_eV"], row["f0"], [process], row)[
            "channel_power"
        ][0]
        self.assertGreater(sdcs, 0.0)
        self.assertEqual(threshold_loss, 0.0)

    def test_quadrature_matches_maxwellian_mean_energy(self):
        temperature = 2.4
        edges = np.linspace(0.0, 30.0 * temperature, 20001)
        energy = (edges[:-1] + edges[1:]) / 2
        f0 = np.exp(-energy / temperature)
        result = _quadrature(
            energy,
            f0,
            [],
            {
                "energy_eV": energy,
                "energy_edges_eV": edges,
                "gas_temperature_K": 300.0,
                "target_fractions": np.array([]),
                "product_fractions": np.array([]),
            },
        )
        self.assertAlmostEqual(result["mean_energy_eV"], 1.5 * temperature, delta=1e-4)


if __name__ == "__main__":
    unittest.main()
