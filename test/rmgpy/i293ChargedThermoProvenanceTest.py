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

"""Regression tests for charged-species thermo provenance in plasma reactors."""

import copy
import math
import os
import pickle
import re

import pytest
import yaml

import rmgpy.data.rmg as rmg_data_module
from rmgpy import constants
from rmgpy import settings
from rmgpy.chemkin import load_chemkin_file, save_chemkin_file, save_species_dictionary
from rmgpy.data.base import Entry
from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.thermo import ThermoLibrary
from rmgpy.exceptions import PlasmaStateError
from rmgpy.molecule import Molecule
from rmgpy.rmg.model import ReactionModel
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial, ThermoData
from rmgpy.thermo.thermoengine import generate_thermo_data, process_thermo_data
from rmgpy.yaml_cantera2 import save_cantera_model


DATABASE_PATH = os.path.join(settings["database.directory"], "thermo")
ARGON_ION_ADJACENCY = "multiplicity 2\n1 Ar u1 p3 c+1"
ARGON_DIMER_ION_ADJACENCY = (
    "multiplicity 2\n"
    "1 Ar u0 p3 c+1 {2,S}\n"
    "2 Ar u1 p3 c0 {1,S}"
)
THERMO_SOURCE_ASSERTION = "caller-asserted, not verified"

pytestmark = pytest.mark.database


@pytest.fixture(scope="module")
def thermo_database():
    previous = rmg_data_module.database
    rmg_data_module.database = None
    try:
        database = RMGDatabase()
        database.load_thermo(
            DATABASE_PATH,
            thermo_libraries=["primaryThermoLibrary", "PlasmaThermo"],
            depository=False,
        )
        yield database.thermo
    finally:
        rmg_data_module.database = previous


def _thermo_data(h298_kj=0.0):
    """Return a real thermo object, never a comment-only stand-in."""
    return ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
        Cpdata=([30.0] * 7, "J/(mol*K)"),
        H298=(h298_kj, "kJ/mol"),
        S298=(180.0, "J/(mol*K)"),
        Cp0=(4.0 * 8.314462618, "J/(mol*K)"),
        CpInf=(7.0 * 8.314462618, "J/(mol*K)"),
    )


def _thermo_data_with_hidden_cp_error():
    """Return a pair that only differs at a missed ThermoData breakpoint."""
    reference = _thermo_data()
    provisional = ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
        Cpdata=([30.0, 31.0, 30.0, 30.0, 30.0, 30.0, 30.0], "J/(mol*K)"),
        H298=(0.0, "J/mol"),
        S298=(180.0, "J/(mol*K)"),
        Cp0=(4.0 * 8.314462618, "J/(mol*K)"),
        CpInf=(7.0 * 8.314462618, "J/(mol*K)"),
    )
    old_reference_temperature = 1150.0
    enthalpy_offset = (
        provisional.get_enthalpy(old_reference_temperature)
        - reference.get_enthalpy(old_reference_temperature)
    )
    entropy_offset = (
        provisional.get_entropy(old_reference_temperature)
        - reference.get_entropy(old_reference_temperature)
    )
    perturbed = ThermoData(
        Tdata=provisional.Tdata,
        Cpdata=provisional.Cpdata,
        H298=(-enthalpy_offset, "J/mol"),
        S298=(180.0 - entropy_offset, "J/(mol*K)"),
        Cp0=provisional.Cp0,
        CpInf=provisional.CpInf,
    )
    return reference, perturbed


def _mixed_form_thermo_with_hidden_cp_error():
    """Return ThermoData and NASA that differ only between table nodes."""
    gas_constant = constants.R
    reference = ThermoData(
        Tdata=([300.0, 1500.0], "K"),
        Cpdata=([3.5 * gas_constant] * 2, "J/(mol*K)"),
        H298=(0.0, "J/mol"),
        S298=(180.0, "J/(mol*K)"),
        Cp0=(3.5 * gas_constant, "J/(mol*K)"),
        CpInf=(3.5 * gas_constant, "J/(mol*K)"),
        Tmin=(300.0, "K"),
        Tmax=(1500.0, "K"),
    )
    reference_temperature = 900.0
    epsilon = -0.12 / 360000.0
    a0 = 3.5 + epsilon * 450000.0
    a1 = -1800.0 * epsilon
    a2 = epsilon
    a5 = reference_temperature * (
        reference.get_enthalpy(reference_temperature)
        / (gas_constant * reference_temperature)
        - a0
        - a1 * reference_temperature / 2.0
        - a2 * reference_temperature ** 2 / 3.0
    )
    a6 = (
        reference.get_entropy(reference_temperature) / gas_constant
        - a0 * math.log(reference_temperature)
        - a1 * reference_temperature
        - a2 * reference_temperature ** 2 / 2.0
    )
    actual = NASA(
        polynomials=[NASAPolynomial(
            coeffs=[a0, a1, a2, 0.0, 0.0, a5, a6],
            Tmin=(300.0, "K"),
            Tmax=(1500.0, "K"),
        )],
        Tmin=(300.0, "K"),
        Tmax=(1500.0, "K"),
    )
    return reference, actual


def _high_temperature_nasa(tmin=2500.0, tmax=3500.0):
    polynomial = NASAPolynomial(
        coeffs=[3.5, 1.0e-4, -2.0e-8, 3.0e-12, -1.0e-16, 1200.0, 4.0],
        Tmin=(tmin, "K"),
        Tmax=(tmax, "K"),
    )
    return NASA(polynomials=[polynomial], Tmin=(tmin, "K"), Tmax=(tmax, "K"))


def _segmented_nasa(narrow_cp=3.5, include_narrow=True):
    segments = [
        (300.0, 1000.0, 3.5),
        (1000.0, 1010.0, narrow_cp),
        (1010.0, 3000.0, 3.5),
    ]
    if not include_narrow:
        del segments[1]
    polynomials = [
        NASAPolynomial(
            coeffs=[cp, 0.0, 0.0, 0.0, 0.0, 1200.0, 4.0],
            Tmin=(tmin, "K"),
            Tmax=(tmax, "K"),
        )
        for tmin, tmax, cp in segments
    ]
    return NASA(polynomials=polynomials, Tmin=(300.0, "K"), Tmax=(3000.0, "K"))


def _cp_nullspace_nasa(perturbed=False):
    """Return the Round 69 NASA-7 counterexample over one polynomial piece."""
    coefficients = [3.5, 0.0, 0.0, 0.0, 0.0, 1200.0, 4.0]
    if perturbed:
        delta = [-0.00374034375, 7.711875e-6, -4.95e-9, 1.0e-12, 0.0]
        temperature = 1650.0
        delta_h = sum(
            coefficient * temperature ** (index + 1) / (index + 1)
            for index, coefficient in enumerate(delta)
        )
        delta_s = (
            delta[0] * math.log(temperature)
            + delta[1] * temperature
            + delta[2] * temperature ** 2 / 2.0
            + delta[3] * temperature ** 3 / 3.0
            + delta[4] * temperature ** 4 / 4.0
        )
        for index, coefficient in enumerate(delta):
            coefficients[index] += coefficient
        coefficients[5] -= delta_h
        coefficients[6] -= delta_s
    polynomial = NASAPolynomial(
        coeffs=coefficients,
        Tmin=(300.0, "K"),
        Tmax=(3000.0, "K"),
    )
    return NASA(polynomials=[polynomial], Tmin=(300.0, "K"), Tmax=(3000.0, "K"))


class _CountingEntries(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.values_calls = 0

    def values(self):
        self.values_calls += 1
        return super().values()


class _NeutralStructureTrap:
    """A neutral library structure that must be skipped before formula work."""

    def get_net_charge(self):
        return 0

    def get_formula(self):
        raise AssertionError("neutral library entries must not be formula-bucketed")


class _MalformedThermo:
    """An entry payload without the heat-capacity model contract."""

    pass


def _electron():
    return Species(label="e-").from_adjacency_list("1 e u1 p0 c-1")


def _initialize(ion, assertions=None, edge_species=None):
    electron = _electron()
    initial = {ion: 0.5, electron: 0.5}
    reactor = PlasmaReactor(
        (300, "K"),
        (1, "bar"),
        initial,
        (300, "K"),
        thermo_source_assertions=assertions,
    )
    reactor.initialize_model([ion, electron], [], edge_species or [], [])
    return reactor


def _library_argon_ion(thermo_database):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    library_thermo = thermo_database.get_thermo_data(ion)
    ion.thermo = process_thermo_data(ion, library_thermo)
    return ion


def _proton():
    return Species(label="H+").from_adjacency_list("1 H u0 p0 c+1")


def _add_library(thermo_database, library):
    thermo_database.libraries[library.label] = library
    thermo_database.library_order.insert(0, library.label)


def _remove_library(thermo_database, library):
    thermo_database.library_order.remove(library.label)
    del thermo_database.libraries[library.label]


def _hbi_estimated_ammonia_ion(thermo_database):
    """Build a genuine library-seeded HBI estimate for NH3+ in memory."""
    ammonium = Molecule().from_adjacency_list(
        "1 N u0 p0 c+1 {2,S} {3,S} {4,S} {5,S}\n"
        "2 H u0 p0 c0 {1,S}\n"
        "3 H u0 p0 c0 {1,S}\n"
        "4 H u0 p0 c0 {1,S}\n"
        "5 H u0 p0 c0 {1,S}"
    )
    library = ThermoLibrary(label="ChargedHBITest", thermo_convention="ion")
    library.entries["NH4+"] = Entry(
        index=1,
        label="NH4+",
        item=ammonium,
        data=_thermo_data(600.0),
    )
    thermo_database.libraries[library.label] = library
    thermo_database.library_order.insert(0, library.label)
    ion = Species(label="NH3+").from_adjacency_list(
        "multiplicity 2\n"
        "1 N u1 p0 c+1 {2,S} {3,S} {4,S}\n"
        "2 H u0 p0 c0 {1,S}\n"
        "3 H u0 p0 c0 {1,S}\n"
        "4 H u0 p0 c0 {1,S}"
    )
    try:
        ion.thermo = thermo_database.get_thermo_data(ion)
    finally:
        thermo_database.library_order.remove(library.label)
        del thermo_database.libraries[library.label]
    return ion


def _from_cantera_species(cantera_species, molecule):
    coefficients = cantera_species.thermo.coeffs
    midpoint = coefficients[0]
    thermo = NASA(
        polynomials=[
            NASAPolynomial(
                coeffs=coefficients[8:15],
                Tmin=(cantera_species.thermo.min_temp, "K"),
                Tmax=(midpoint, "K"),
            ),
            NASAPolynomial(
                coeffs=coefficients[1:8],
                Tmin=(midpoint, "K"),
                Tmax=(cantera_species.thermo.max_temp, "K"),
            ),
        ],
        Tmin=(cantera_species.thermo.min_temp, "K"),
        Tmax=(cantera_species.thermo.max_temp, "K"),
    )
    return Species(label=cantera_species.name, molecule=[copy.deepcopy(molecule)], thermo=thermo)


def test_library_seeded_hbi_with_library_comment_is_refused(thermo_database):
    ion = _hbi_estimated_ammonia_ion(thermo_database)
    assert ion.thermo.comment.startswith("Thermo library: ChargedHBITest + radical(")
    with pytest.raises(PlasmaStateError, match="could not be value-matched"):
        _initialize(ion)


def test_library_ion_is_accepted_after_pickle_restart(thermo_database):
    restored = pickle.loads(pickle.dumps(_library_argon_ion(thermo_database)))
    _initialize(restored)


def test_electrochemical_convention_proton_is_refused_by_name(thermo_database):
    proton = _proton()
    library = ThermoLibrary(label="ElectrochemicalFixture", thermo_convention="ion")
    library.thermo_convention = "electrochemical"
    library.entries["proton"] = Entry(
        index=1, label="proton", item=proton.molecule[0], data=_thermo_data(0.0)
    )
    _add_library(thermo_database, library)
    proton.thermo = process_thermo_data(proton, copy.deepcopy(library.entries["proton"].data))
    try:
        with pytest.raises(PlasmaStateError) as error:
            _initialize(proton)
        message = str(error.value)
        assert "H+" in message
        assert "ElectrochemicalFixture/proton" in message
        assert "electrochemical" in message
    finally:
        _remove_library(thermo_database, library)


def test_ion_convention_proton_is_accepted(thermo_database):
    proton = _proton()
    library = ThermoLibrary(label="IonConventionFixture", thermo_convention="ion")
    library.thermo_convention = "ion"
    library.entries["proton"] = Entry(
        index=1, label="proton", item=proton.molecule[0], data=_thermo_data(1530.0)
    )
    _add_library(thermo_database, library)
    proton.thermo = process_thermo_data(proton, copy.deepcopy(library.entries["proton"].data))
    try:
        _initialize(proton)
    finally:
        _remove_library(thermo_database, library)


def test_electron_from_electrocat_thermo_is_accepted(thermo_database):
    electron = _electron()
    library = ThermoLibrary(label="electrocatThermo", thermo_convention="ion")
    library.thermo_convention = "electrochemical"
    library.entries["electron"] = Entry(
        index=1, label="electron", item=electron.molecule[0], data=_thermo_data()
    )
    _add_library(thermo_database, library)
    electron.thermo = process_thermo_data(
        electron, copy.deepcopy(library.entries["electron"].data)
    )
    neutral = Species(label="N2", smiles="N#N")
    neutral.thermo = thermo_database.get_thermo_data(neutral)
    reactor = PlasmaReactor(
        (300, "K"), (1, "bar"), {neutral: 1.0, electron: 1.0e-12}, (300, "K")
    )
    try:
        reactor.initialize_model([neutral, electron], [], [], [])
    finally:
        _remove_library(thermo_database, library)


def test_library_ion_is_accepted_after_nonverbose_chemkin_reload(thermo_database, tmp_path):
    ion = _library_argon_ion(thermo_database)
    chemkin = tmp_path / "chem.inp"
    dictionary = tmp_path / "species_dictionary.txt"
    save_species_dictionary(str(dictionary), [ion])
    save_chemkin_file(str(chemkin), [ion], [], verbose=False)

    loaded_species, _ = load_chemkin_file(str(chemkin), str(dictionary))
    loaded = loaded_species[0]
    assert loaded.thermo.comment == ""
    _initialize(loaded)


def test_library_ion_is_accepted_after_notes_stripped_cantera_reload(
    thermo_database, tmp_path
):
    import cantera as ct

    ion = _library_argon_ion(thermo_database)
    yaml_path = tmp_path / "chem.yaml"
    stripped_path = tmp_path / "chem-no-notes.yaml"
    save_cantera_model(ReactionModel(species=[ion], reactions=[]), str(yaml_path))
    document = yaml.safe_load(yaml_path.read_text())
    for species in document["species"]:
        species.pop("note", None)
        species["thermo"].pop("note", None)
    document["phases"][0]["transport"] = "none"
    stripped_path.write_text(yaml.safe_dump(document, sort_keys=False))

    gas = ct.Solution(str(stripped_path))
    cantera_species = gas.species(0)
    assert "note" not in cantera_species.input_data
    assert "note" not in cantera_species.input_data["thermo"]
    loaded = _from_cantera_species(cantera_species, ion.molecule[0])
    _initialize(loaded)


def test_library_entry_with_prior_comment_is_accepted(thermo_database):
    entry = thermo_database.libraries["PlasmaThermo"].entries["[Arp]"]
    original_comment = entry.data.comment
    entry.data.comment = "Curated JANAF ion data. "
    try:
        ion = _library_argon_ion(thermo_database)
    finally:
        entry.data.comment = original_comment
    assert ion.thermo.comment.startswith("Curated JANAF ion data.")
    _initialize(ion)


def test_library_thermodata_without_cp_limits_is_accepted(thermo_database):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    data = _thermo_data(1500.0)
    data.Cp0 = None
    data.CpInf = None
    library = ThermoLibrary(label="MissingCpLimits", thermo_convention="ion")
    library.entries["Ar+"] = Entry(
        index=1, label="Ar+", item=ion.molecule[0], data=data
    )
    thermo_database.libraries[library.label] = library
    thermo_database.library_order.insert(0, library.label)
    try:
        ion.thermo = process_thermo_data(ion, thermo_database.get_thermo_data(ion))
        assert library.entries["Ar+"].data.Cp0 is None
        assert library.entries["Ar+"].data.CpInf is None
        _initialize(ion)
    finally:
        thermo_database.library_order.remove(library.label)
        del thermo_database.libraries[library.label]


def test_real_library_thermodata_processed_by_normal_engine_is_accepted(
    thermo_database,
):
    ion = Species(label="Ar2p").from_adjacency_list(
        ARGON_DIMER_ION_ADJACENCY
    )
    library_thermo = thermo_database.get_thermo_data(ion)
    assert isinstance(library_thermo, ThermoData)
    assert library_thermo.label == "[Ar2p]"
    ion.thermo = generate_thermo_data(ion)
    assert isinstance(ion.thermo, NASA)
    _initialize(ion)


def test_thermodata_cp_error_at_breakpoint_is_refused(thermo_database):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    reference, perturbed = _thermo_data_with_hidden_cp_error()
    library = ThermoLibrary(label="ThermoDataBreakpointIon", thermo_convention="ion")
    library.entries["Ar+-thermodata"] = Entry(
        index=1,
        label="Ar+-thermodata",
        item=ion.molecule[0],
        data=reference,
    )
    thermo_database.libraries[library.label] = library
    thermo_database.library_order.insert(0, library.label)
    ion.thermo = perturbed
    cp_difference = (
        ion.thermo.get_heat_capacity(400.0)
        - reference.get_heat_capacity(400.0)
    )
    assert cp_difference == pytest.approx(1.0)
    assert ion.thermo.get_enthalpy(1150.0) == pytest.approx(
        reference.get_enthalpy(1150.0), rel=0.0, abs=1.0e-9
    )
    assert ion.thermo.get_entropy(1150.0) == pytest.approx(
        reference.get_entropy(1150.0), rel=0.0, abs=1.0e-10
    )
    try:
        with pytest.raises(
            PlasmaStateError, match="could not be value-matched"
        ):
            _initialize(ion)
    finally:
        thermo_database.library_order.remove(library.label)
        del thermo_database.libraries[library.label]


def test_mixed_nasa_thermodata_cp_error_between_nodes_is_refused(
    thermo_database,
):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    reference, actual = _mixed_form_thermo_with_hidden_cp_error()
    library = ThermoLibrary(label="MixedFormIon", thermo_convention="ion")
    library.entries["Ar+-thermodata"] = Entry(
        index=1,
        label="Ar+-thermodata",
        item=ion.molecule[0],
        data=reference,
    )
    thermo_database.libraries[library.label] = library
    thermo_database.library_order.insert(0, library.label)
    ion.thermo = actual
    assert actual.get_heat_capacity(300.0) == pytest.approx(
        reference.get_heat_capacity(300.0)
    )
    assert actual.get_heat_capacity(1500.0) == pytest.approx(
        reference.get_heat_capacity(1500.0)
    )
    assert actual.get_heat_capacity(900.0) - reference.get_heat_capacity(
        900.0
    ) == pytest.approx(0.99773664)
    assert actual.get_enthalpy(900.0) == pytest.approx(
        reference.get_enthalpy(900.0), abs=1.0e-9
    )
    assert actual.get_entropy(900.0) == pytest.approx(
        reference.get_entropy(900.0), abs=1.0e-10
    )
    try:
        with pytest.raises(
            PlasmaStateError, match="could not be value-matched"
        ):
            _initialize(ion)
    finally:
        thermo_database.library_order.remove(library.label)
        del thermo_database.libraries[library.label]


def test_wilhoit_match_diagnostic_states_that_comparison_is_sampled(
    thermo_database,
):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    reference = _thermo_data().to_wilhoit(B=1000.0)
    library = ThermoLibrary(label="WilhoitIon", thermo_convention="ion")
    library.entries["Ar+-wilhoit"] = Entry(
        index=1,
        label="Ar+-wilhoit",
        item=ion.molecule[0],
        data=reference,
    )
    thermo_database.libraries[library.label] = library
    thermo_database.library_order.insert(0, library.label)
    ion.thermo = copy.deepcopy(reference)
    try:
        reactor = _initialize(ion)
        diagnostic = reactor.thermo_provenance_diagnostics["Ar+"]
        assert "value-matched to library WilhoitIon/Ar+-wilhoit" in diagnostic
        assert "sampled Wilhoit comparison" in diagnostic
        assert "20 log-spaced Cp points" in diagnostic
    finally:
        thermo_database.library_order.remove(library.label)
        del thermo_database.libraries[library.label]


def test_exact_library_nasa_valid_only_above_2000_k_is_accepted(thermo_database):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    library = ThermoLibrary(label="HighTemperatureIon", thermo_convention="ion")
    library.entries["Ar+-hot"] = Entry(
        index=1,
        label="Ar+-hot",
        item=ion.molecule[0],
        data=_high_temperature_nasa(),
    )
    thermo_database.libraries[library.label] = library
    thermo_database.library_order.insert(0, library.label)
    ion.thermo = copy.deepcopy(library.entries["Ar+-hot"].data)
    try:
        _initialize(ion)
    finally:
        thermo_database.library_order.remove(library.label)
        del thermo_database.libraries[library.label]


def test_models_with_no_overlapping_valid_temperature_are_not_matched(thermo_database):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    library = ThermoLibrary(label="DisjointTemperatureIon", thermo_convention="ion")
    library.entries["Ar+-cold"] = Entry(
        index=1,
        label="Ar+-cold",
        item=ion.molecule[0],
        data=_high_temperature_nasa(2500.0, 3500.0),
    )
    thermo_database.libraries[library.label] = library
    thermo_database.library_order.insert(0, library.label)
    ion.thermo = _high_temperature_nasa(4000.0, 5000.0)
    try:
        with pytest.raises(PlasmaStateError, match="could not be value-matched"):
            _initialize(ion)
    finally:
        thermo_database.library_order.remove(library.label)
        del thermo_database.libraries[library.label]


def test_narrow_nasa_segment_mismatch_is_refused(thermo_database):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    library = ThermoLibrary(label="NarrowSegmentIon", thermo_convention="ion")
    library.entries["Ar+-segmented"] = Entry(
        index=1,
        label="Ar+-segmented",
        item=ion.molecule[0],
        data=_segmented_nasa(),
    )
    thermo_database.libraries[library.label] = library
    thermo_database.library_order.insert(0, library.label)
    ion.thermo = _segmented_nasa(narrow_cp=4.5)
    try:
        with pytest.raises(PlasmaStateError, match="could not be value-matched"):
            _initialize(ion)
    finally:
        thermo_database.library_order.remove(library.label)
        del thermo_database.libraries[library.label]


def test_degree_four_cp_nullspace_mismatch_is_refused(thermo_database):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    library = ThermoLibrary(label="CpNullspaceIon", thermo_convention="ion")
    reference = _cp_nullspace_nasa()
    library.entries["Ar+-counterexample"] = Entry(
        index=1,
        label="Ar+-counterexample",
        item=ion.molecule[0],
        data=reference,
    )
    thermo_database.libraries[library.label] = library
    thermo_database.library_order.insert(0, library.label)
    ion.thermo = _cp_nullspace_nasa(perturbed=True)
    for temperature in (975.0, 1650.0, 2325.0):
        assert ion.thermo.get_heat_capacity(temperature) == pytest.approx(
            reference.get_heat_capacity(temperature), rel=0.0, abs=1.0e-10
        )
    assert ion.thermo.get_enthalpy(1650.0) == pytest.approx(
        reference.get_enthalpy(1650.0), rel=0.0, abs=1.0e-9
    )
    assert ion.thermo.get_entropy(1650.0) == pytest.approx(
        reference.get_entropy(1650.0), rel=0.0, abs=1.0e-10
    )
    cp_difference = abs(
        ion.thermo.get_heat_capacity(500.0)
        - reference.get_heat_capacity(500.0)
    )
    assert cp_difference == pytest.approx(0.00829, rel=1.0e-3)
    try:
        with pytest.raises(PlasmaStateError, match="could not be value-matched"):
            _initialize(ion)
    finally:
        thermo_database.library_order.remove(library.label)
        del thermo_database.libraries[library.label]


def test_nasa_coverage_gap_is_a_contextual_plasma_refusal(thermo_database):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    library = ThermoLibrary(label="GappedIon", thermo_convention="ion")
    library.entries["Ar+-gap"] = Entry(
        index=1,
        label="Ar+-gap",
        item=ion.molecule[0],
        data=_segmented_nasa(include_narrow=False),
    )
    thermo_database.libraries[library.label] = library
    thermo_database.library_order.insert(0, library.label)
    ion.thermo = _segmented_nasa()
    try:
        with pytest.raises(PlasmaStateError) as error:
            _initialize(ion)
        message = str(error.value)
        assert "Ar+" in message
        assert "GappedIon/Ar+-gap" in message
        assert "coverage gap" in message
        assert "1000" in message and "1010" in message
    finally:
        thermo_database.library_order.remove(library.label)
        del thermo_database.libraries[library.label]


def test_nasa_gap_outside_candidate_overlap_is_a_contextual_refusal(thermo_database):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    library = ThermoLibrary(label="ShortRangeIon", thermo_convention="ion")
    short_polynomial = NASAPolynomial(
        coeffs=[3.5, 0.0, 0.0, 0.0, 0.0, 1200.0, 4.0],
        Tmin=(300.0, "K"),
        Tmax=(900.0, "K"),
    )
    library.entries["Ar+-short"] = Entry(
        index=1,
        label="Ar+-short",
        item=ion.molecule[0],
        data=NASA(
            polynomials=[short_polynomial],
            Tmin=(300.0, "K"),
            Tmax=(900.0, "K"),
        ),
    )
    thermo_database.libraries[library.label] = library
    thermo_database.library_order.insert(0, library.label)
    ion.thermo = _segmented_nasa(include_narrow=False)
    try:
        with pytest.raises(PlasmaStateError) as error:
            _initialize(ion)
        message = str(error.value)
        assert "Ar+" in message
        assert "ShortRangeIon/Ar+-short" in message
        assert "species thermo" in message
        assert "coverage gap" in message
        assert "1000" in message and "1010" in message
    finally:
        thermo_database.library_order.remove(library.label)
        del thermo_database.libraries[library.label]


def test_gapped_charged_edge_is_refused_before_and_after_promotion(thermo_database):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    library = ThermoLibrary(label="EdgeGapIon", thermo_convention="ion")
    library.entries["Ar+-short"] = Entry(
        index=1,
        label="Ar+-short",
        item=ion.molecule[0],
        data=NASA(
            polynomials=[NASAPolynomial(
                coeffs=[3.5, 0.0, 0.0, 0.0, 0.0, 1200.0, 4.0],
                Tmin=(300.0, "K"),
                Tmax=(900.0, "K"),
            )],
            Tmin=(300.0, "K"),
            Tmax=(900.0, "K"),
        ),
    )
    thermo_database.libraries[library.label] = library
    thermo_database.library_order.insert(0, library.label)
    ion.thermo = _segmented_nasa(include_narrow=False)
    neutral = Species(label="Ar").from_adjacency_list("1 Ar u0 p4 c0")
    neutral.thermo = thermo_database.get_thermo_data(neutral)
    electron = _electron()
    reactor = PlasmaReactor(
        (300, "K"),
        (1, "bar"),
        {neutral: 1.0 - 1.0e-12, electron: 1.0e-12},
        (300, "K"),
    )
    try:
        with pytest.raises(PlasmaStateError) as edge_error:
            reactor.initialize_model([neutral, electron], [], [ion], [])
        message = str(edge_error.value)
        assert "edge species" in message and "Ar+" in message
        assert "EdgeGapIon/Ar+-short" in message
        assert "coverage gap" in message
        assert "1000" in message and "1010" in message

        with pytest.raises(PlasmaStateError) as error:
            reactor.initialize_model([neutral, electron, ion], [], [], [])
        message = str(error.value)
        assert "EdgeGapIon/Ar+-short" in message
        assert "coverage gap" in message
        assert "1000" in message and "1010" in message
    finally:
        thermo_database.library_order.remove(library.label)
        del thermo_database.libraries[library.label]


def test_gapped_earlier_candidate_is_recorded_and_later_match_is_used(
    thermo_database,
):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    earlier = ThermoLibrary(label="GappedEarlier", thermo_convention="ion")
    earlier.entries["Ar+-gap"] = Entry(
        index=1,
        label="Ar+-gap",
        item=ion.molecule[0],
        data=_segmented_nasa(include_narrow=False),
    )
    later = ThermoLibrary(label="ValidLater", thermo_convention="ion")
    later.entries["Ar+-valid"] = Entry(
        index=2,
        label="Ar+-valid",
        item=ion.molecule[0],
        data=_segmented_nasa(),
    )
    for library in (earlier, later):
        thermo_database.libraries[library.label] = library
    labels = [earlier.label, later.label]
    thermo_database.library_order[0:0] = labels
    ion.thermo = _segmented_nasa()
    try:
        reactor = _initialize(ion)
        diagnostic = reactor.thermo_provenance_diagnostics["Ar+"]
        assert "value-matched to library ValidLater/Ar+-valid" in diagnostic
        assert "GappedEarlier/Ar+-gap" in diagnostic
        assert "coverage gap" in diagnostic
    finally:
        del thermo_database.library_order[:len(labels)]
        for label in labels:
            del thermo_database.libraries[label]


def test_malformed_earlier_candidate_is_recorded_and_later_match_is_used(
    thermo_database,
):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    earlier = ThermoLibrary(label="MalformedEarlier", thermo_convention="ion")
    earlier.entries["Ar+-broken"] = Entry(
        index=1,
        label="Ar+-broken",
        item=ion.molecule[0],
        data=_MalformedThermo(),
    )
    later = ThermoLibrary(label="ValidAfterMalformed", thermo_convention="ion")
    later.entries["Ar+-valid"] = Entry(
        index=2,
        label="Ar+-valid",
        item=ion.molecule[0],
        data=_segmented_nasa(),
    )
    for library in (earlier, later):
        thermo_database.libraries[library.label] = library
    labels = [earlier.label, later.label]
    thermo_database.library_order[0:0] = labels
    ion.thermo = _segmented_nasa()
    try:
        reactor = _initialize(ion)
        diagnostic = reactor.thermo_provenance_diagnostics["Ar+"]
        assert (
            "value-matched to library ValidAfterMalformed/Ar+-valid"
            in diagnostic
        )
        assert "MalformedEarlier/Ar+-broken" in diagnostic
        assert "malformed" in diagnostic
    finally:
        del thermo_database.library_order[:len(labels)]
        for label in labels:
            del thermo_database.libraries[label]


def test_liquid_libraries_are_excluded_and_first_gas_match_is_reported(thermo_database):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    data = _thermo_data(1500.0)
    libraries = [
        ThermoLibrary(label="LiquidMatch", thermo_convention="ion"),
        ThermoLibrary(label="GasFirst", thermo_convention="ion"),
        ThermoLibrary(label="GasSecond", thermo_convention="ion"),
    ]
    libraries[0].solvent = "water"
    for index, library in enumerate(libraries, start=1):
        library.entries["Ar+"] = Entry(
            index=index,
            label="Ar+",
            item=ion.molecule[0],
            data=copy.deepcopy(data),
        )
        thermo_database.libraries[library.label] = library
    labels = [library.label for library in libraries]
    thermo_database.library_order[0:0] = labels
    ion.thermo = process_thermo_data(ion, copy.deepcopy(data))
    try:
        reactor = _initialize(ion)
        assert reactor.thermo_provenance_diagnostics["Ar+"] == (
            "value-matched to library GasFirst/Ar+ (first match in library_order)"
        )
    finally:
        del thermo_database.library_order[:len(labels)]
        for label in labels:
            del thermo_database.libraries[label]


def test_library_index_is_built_once_per_initialize_model_call(thermo_database):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    data = _thermo_data(1500.0)
    library = ThermoLibrary(label="CountedLibrary", thermo_convention="ion")
    library.entries = _CountingEntries({
        "Ar+": Entry(index=1, label="Ar+", item=ion.molecule[0], data=data),
    })
    thermo_database.libraries[library.label] = library
    thermo_database.library_order.insert(0, library.label)
    ion.thermo = process_thermo_data(ion, copy.deepcopy(data))
    electron = _electron()
    core = [ion, electron]
    reactor = PlasmaReactor(
        (300, "K"), (1, "bar"), {ion: 0.5, electron: 0.5}, (300, "K")
    )
    try:
        reactor.initialize_model(core, [], [], [])
        assert library.entries.values_calls == 1
        reactor.initialize_model(core, [], [], [])
        assert library.entries.values_calls == 2
    finally:
        thermo_database.library_order.remove(library.label)
        del thermo_database.libraries[library.label]


def test_library_index_skips_neutral_entries_before_formula_bucketing(thermo_database):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    data = _thermo_data(1500.0)
    library = ThermoLibrary(label="ChargedEntriesOnly", thermo_convention="ion")
    library.entries["neutral"] = Entry(
        index=1,
        label="neutral",
        item=_NeutralStructureTrap(),
        data=_thermo_data(),
    )
    library.entries["Ar+"] = Entry(
        index=2,
        label="Ar+",
        item=ion.molecule[0],
        data=data,
    )
    thermo_database.libraries[library.label] = library
    thermo_database.library_order.insert(0, library.label)
    ion.thermo = process_thermo_data(ion, copy.deepcopy(data))
    try:
        _initialize(ion)
    finally:
        thermo_database.library_order.remove(library.label)
        del thermo_database.libraries[library.label]


def test_in_place_species_thermo_mutation_changes_next_verdict(thermo_database):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    data = _thermo_data(1500.0)
    library = ThermoLibrary(label="MutableSpeciesThermo", thermo_convention="ion")
    library.entries["Ar+"] = Entry(
        index=1, label="Ar+", item=ion.molecule[0], data=data
    )
    thermo_database.libraries[library.label] = library
    thermo_database.library_order.insert(0, library.label)
    ion.thermo = process_thermo_data(ion, copy.deepcopy(data))
    electron = _electron()
    core = [ion, electron]
    reactor = PlasmaReactor(
        (300, "K"), (1, "bar"), {ion: 0.5, electron: 0.5}, (300, "K")
    )
    try:
        reactor.initialize_model(core, [], [], [])
        polynomial = ion.thermo.polynomials[0]
        coefficients = polynomial.coeffs
        coefficients[0] += 1.0
        polynomial.coeffs = coefficients
        with pytest.raises(PlasmaStateError, match="could not be value-matched"):
            reactor.initialize_model(core, [], [], [])
    finally:
        thermo_database.library_order.remove(library.label)
        del thermo_database.libraries[library.label]


def test_same_length_library_entry_replacement_changes_next_verdict(thermo_database):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    data = _thermo_data(1500.0)
    library = ThermoLibrary(label="ReplaceableLibraryEntry", thermo_convention="ion")
    library.entries["Ar+"] = Entry(
        index=1, label="Ar+", item=ion.molecule[0], data=data
    )
    thermo_database.libraries[library.label] = library
    thermo_database.library_order.insert(0, library.label)
    ion.thermo = process_thermo_data(ion, copy.deepcopy(data))
    electron = _electron()
    core = [ion, electron]
    reactor = PlasmaReactor(
        (300, "K"), (1, "bar"), {ion: 0.5, electron: 0.5}, (300, "K")
    )
    try:
        reactor.initialize_model(core, [], [], [])
        library.entries["Ar+"] = Entry(
            index=2,
            label="Ar+",
            item=ion.molecule[0],
            data=_thermo_data(1700.0),
        )
        assert len(library.entries) == 1
        with pytest.raises(PlasmaStateError, match="could not be value-matched"):
            reactor.initialize_model(core, [], [], [])
    finally:
        thermo_database.library_order.remove(library.label)
        del thermo_database.libraries[library.label]


@pytest.mark.parametrize(
    "label,smiles",
    [("N2+", "[N+]#N"), ("NO+", "[N+]=O"), ("O2-", "[O-][O]")],
)
def test_gav_ions_are_refused_in_core(thermo_database, label, smiles):
    ion = Species(label=label).from_smiles(smiles)
    ion.thermo = thermo_database.get_thermo_data(ion)
    assert "group additivity" in ion.thermo.comment
    with pytest.raises(PlasmaStateError, match=re.escape(label)):
        _initialize(ion)


def test_electron_is_exempt_without_thermo(thermo_database):
    electron = _electron()
    neutral = Species(label="N2", smiles="N#N")
    neutral.thermo = thermo_database.get_thermo_data(neutral)
    reactor = PlasmaReactor((300, "K"), (1, "bar"), {neutral: 1.0, electron: 1.0e-12}, (300, "K"))
    reactor.initialize_model([neutral, electron], [], [], [])


def test_gav_neutral_is_accepted(thermo_database):
    neutral = Species(label="CH3F", smiles="CF")
    neutral.thermo = thermo_database.get_thermo_data(neutral)
    assert "group additivity" in neutral.thermo.comment
    electron = _electron()
    reactor = PlasmaReactor((300, "K"), (1, "bar"), {neutral: 1.0, electron: 1.0e-12}, (300, "K"))
    reactor.initialize_model([neutral, electron], [], [], [])


def test_unmatched_edge_ion_is_refused_before_and_after_promotion(thermo_database):
    ion = Species(label="N2+", smiles="[N+]#N")
    ion.thermo = thermo_database.get_thermo_data(ion)
    neutral = Species(label="N2", smiles="N#N")
    neutral.thermo = thermo_database.get_thermo_data(neutral)
    electron = _electron()
    reactor = PlasmaReactor((300, "K"), (1, "bar"), {neutral: 1.0, electron: 1.0e-12}, (300, "K"))

    for _ in range(2):
        with pytest.raises(PlasmaStateError, match=r"edge species.*N2\+"):
            reactor.initialize_model([neutral, electron], [], [ion], [])
        assert not reactor._plasma_validated
    with pytest.raises(PlasmaStateError, match=r"N2\+"):
        reactor.initialize_model([neutral, electron, ion], [], [], [])


def test_initialize_model_clears_stale_provenance_diagnostics(thermo_database):
    ion = _library_argon_ion(thermo_database)
    neutral = Species(label="Ar").from_adjacency_list("1 Ar u0 p4 c0")
    neutral.thermo = thermo_database.get_thermo_data(neutral)
    electron = _electron()
    core = [neutral, electron]
    reactor = PlasmaReactor(
        (300, "K"),
        (1, "bar"),
        {neutral: 1.0 - 1.0e-12, electron: 1.0e-12},
        (300, "K"),
    )

    reactor.initialize_model(core, [], [ion], [])
    assert reactor.thermo_provenance_diagnostics["Ar+"].startswith(
        "value-matched to library"
    )
    reactor.initialize_model(core, [], [], [])
    assert reactor.thermo_provenance_diagnostics == {}


def test_no_database_refuses_without_caller_assertion(thermo_database):
    ion = Species(label="Ar+", thermo=_thermo_data(1500.0)).from_adjacency_list(ARGON_ION_ADJACENCY)
    saved = rmg_data_module.database
    rmg_data_module.database = None
    try:
        with pytest.raises(PlasmaStateError, match="no thermo database is loaded.*thermo_source_assertions"):
            _initialize(ion)
    finally:
        rmg_data_module.database = saved


def test_no_database_accepts_explicit_caller_assertion_and_records_diagnostic(thermo_database):
    ion = Species(label="Ar+", thermo=_thermo_data(1500.0)).from_adjacency_list(ARGON_ION_ADJACENCY)
    saved = rmg_data_module.database
    rmg_data_module.database = None
    try:
        reactor = _initialize(ion, assertions={"Ar+": "ion"})
    finally:
        rmg_data_module.database = saved
    assert reactor.thermo_provenance_diagnostics["Ar+"] == THERMO_SOURCE_ASSERTION


def test_deck_declaration_can_name_a_promoted_species_and_survives_writing(monkeypatch):
    import rmgpy.rmg.input as input_module

    electron = _electron()
    neutral = Species(label="Ar").from_adjacency_list("1 Ar u0 p4 c0")
    job = type("Job", (), {"reaction_systems": []})()
    monkeypatch.setattr(input_module, "rmg", job)
    monkeypatch.setattr(input_module, "species_dict", {"e-": electron, "Ar": neutral})

    input_module.plasma_reactor(
        temperature=(300, "K"),
        pressure=(1, "bar"),
        electronTemperature=(300, "K"),
        initialMoleFractions={"Ar": 1.0 - 1.0e-12, "e-": 1.0e-12},
        terminationTime=(1, "s"),
        thermoSourceAssertions={"promoted-ion": "ion"},
    )

    reactor = job.reaction_systems[0]
    assert reactor.thermo_source_assertions == {"promoted-ion": "ion"}
    assert "thermoSourceAssertions = {'promoted-ion': 'ion'}" in input_module._format_plasma_wall(reactor)


def test_none_thermo_fails_closed(thermo_database):
    ion = Species(label="Ar+").from_adjacency_list(ARGON_ION_ADJACENCY)
    with pytest.raises(PlasmaStateError, match="has no thermo data"):
        _initialize(ion)
