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

"""File-based regressions for every charged-thermo convention bypass."""

import copy
from pathlib import Path
from types import SimpleNamespace

import pytest

import rmgpy.data.rmg as rmg_data_module
from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary
from rmgpy.exceptions import PlasmaStateError
from rmgpy.kinetics import Arrhenius
from rmgpy.reaction import Reaction
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.species import Species


FIXTURES = Path(__file__).parent / "test_data/i315_thermo_convention"


def _load(name, library=None):
    library = ThermoLibrary(label=Path(name).stem) if library is None else library
    context = ThermoDatabase()
    return library.load(str(FIXTURES / name), context.local_context, context.global_context)


def _species(library, label, identity=None):
    entry = library.entries[label]
    return Species(label=identity or label, molecule=[copy.deepcopy(entry.item)],
                   thermo=copy.deepcopy(entry.data))


def _database(monkeypatch, library):
    thermo = ThermoDatabase()
    thermo.libraries = {library.label: library}
    thermo.library_order = [library.label]
    monkeypatch.setattr(rmg_data_module, "database", SimpleNamespace(thermo=thermo))
    return thermo


def _electron():
    return Species(label="e-").from_adjacency_list("1 e u1 p0 c-1")


def _initialize_proton(library, assertions=None):
    proton = _species(library, "proton", "H+")
    electron = _electron()
    reactor = PlasmaReactor((300, "K"), (1, "bar"), {proton: 0.5, electron: 0.5},
                            (300, "K"), thermo_source_assertions=assertions)
    reactor.initialize_model([proton, electron], [], [], [])
    return reactor


def test_default_constructor_and_undeclared_file_fail_closed(monkeypatch):
    library = _load("undeclared.py")
    assert library.thermo_convention is None
    _database(monkeypatch, library)
    with pytest.raises(PlasmaStateError, match="undeclared"):
        _initialize_proton(library)


def test_reused_library_resets_previous_file_declaration(monkeypatch):
    context = ThermoDatabase()
    library = _load("ion.py")
    assert library.thermo_convention == "ion"
    _database(monkeypatch, library)
    _initialize_proton(library)
    library.load(str(FIXTURES / "undeclared.py"), context.local_context, context.global_context)
    assert library.thermo_convention is None
    with pytest.raises(PlasmaStateError, match="undeclared"):
        _initialize_proton(library)


def test_external_plasma_thermo_impostor_refused(monkeypatch):
    database = ThermoDatabase()
    database.load_libraries(str(FIXTURES), libraries=[str(FIXTURES / "PlasmaThermo.py")])
    monkeypatch.setattr(rmg_data_module, "database", SimpleNamespace(thermo=database))
    with pytest.raises(PlasmaStateError, match="PlasmaThermo/proton.*undeclared"):
        _initialize_proton(database.libraries["PlasmaThermo"])


def test_electrochemical_edge_refused_before_reverse_growth_rate(monkeypatch):
    library = _load("electrochemical.py")
    _database(monkeypatch, library)
    proton = _species(library, "proton", "H+")
    oxide = _species(library, "oxide", "O-")
    hydrogen = _species(library, "hydrogen", "H")
    oxygen = _species(library, "oxygen", "O")
    electron = _electron()
    reaction = Reaction(reactants=[proton, oxide], products=[hydrogen, oxygen],
                        reversible=True, kinetics=Arrhenius(A=(1, "m^3/(mol*s)")))
    reactor = PlasmaReactor((300, "K"), (1, "bar"),
                            {hydrogen: 0.5, oxygen: 0.5, electron: 1e-12}, (300, "K"))
    with pytest.raises(PlasmaStateError, match="edge species.*H\+.*electrochemical"):
        reactor.initialize_model([hydrogen, oxygen, electron], [], [proton, oxide], [reaction])
    assert reactor.num_edge_reactions == -1


@pytest.mark.parametrize("convention", ["ion", "electrochemical"])
def test_save_reload_preserves_file_declaration(convention, tmp_path):
    library = _load(convention + ".py")
    saved = tmp_path / "saved.py"
    library.save(str(saved))
    assert "thermoConvention = " + repr(convention) in saved.read_text()
    context = ThermoDatabase()
    reloaded = ThermoLibrary().load(str(saved), context.local_context, context.global_context)
    assert reloaded.thermo_convention == convention
    assert reloaded.entries["proton"].data.get_enthalpy(298.15) == pytest.approx(
        library.entries["proton"].data.get_enthalpy(298.15))


def test_standalone_label_only_proton_assertion_refused(monkeypatch):
    monkeypatch.setattr(rmg_data_module, "database", None)
    with pytest.raises(PlasmaStateError, match="H\+.*explicitly assert the ion thermo convention"):
        _initialize_proton(_load("electrochemical.py"), assertions=["H+"])


def test_file_loading_matrix_declared_undeclared_reused_and_roundtrip(tmp_path):
    """Must-fix 6: exercise declarations exclusively through library files."""
    context = ThermoDatabase()
    library = _load("ion.py")
    assert library.thermo_convention == "ion"
    library.load(str(FIXTURES / "undeclared.py"), context.local_context, context.global_context)
    assert library.thermo_convention is None
    library.load(str(FIXTURES / "electrochemical.py"), context.local_context, context.global_context)
    saved = tmp_path / "PlasmaThermo.py"
    library.save(str(saved))
    reloaded = ThermoLibrary().load(str(saved), context.local_context, context.global_context)
    assert reloaded.thermo_convention == "electrochemical"


def test_parsed_ion_declaration_accepts_gas_phase_proton(monkeypatch):
    library = _load("ion.py")
    assert library.thermo_convention == "ion"
    assert library.entries["proton"].data.get_enthalpy(298.15) == pytest.approx(1530000)
    _database(monkeypatch, library)
    _initialize_proton(library)


@pytest.mark.parametrize("selected", [None, ["undeclared"]])
def test_discovery_and_selected_loaders_leave_undeclared_untrusted(selected):
    database = ThermoDatabase()
    database.load_libraries(str(FIXTURES), libraries=selected)
    assert database.libraries["undeclared"].thermo_convention is None


def test_reviewed_content_bridge_ignores_name_and_expires_on_declaration(tmp_path):
    context = ThermoDatabase()
    content = (FIXTURES / "reviewed_legacy.py").read_text()
    renamed = tmp_path / "arbitrary.py"
    renamed.write_text(content)
    library = ThermoLibrary(thermo_convention=None).load(
        str(renamed), context.local_context, context.global_context)
    assert library.thermo_convention == "ion"
    renamed.write_text(content + "\nthermoConvention = 'electrochemical'\n")
    library.load(str(renamed), context.local_context, context.global_context)
    assert library.thermo_convention == "electrochemical"
    renamed.write_text(content + "\n")
    library.load(str(renamed), context.local_context, context.global_context)
    assert library.thermo_convention is None


def test_standalone_explicit_ion_assertion_survives_pickle_and_deck_write(monkeypatch):
    import pickle
    import rmgpy.rmg.input as input_module

    monkeypatch.setattr(rmg_data_module, "database", None)
    reactor = _initialize_proton(_load("ion.py"), assertions={"H+": "ion"})
    restored = pickle.loads(pickle.dumps(reactor))
    assert restored.thermo_source_assertions == {"H+": "ion"}
    assert "thermoSourceAssertions = {'H+': 'ion'}" in input_module._format_plasma_wall(restored)


def test_standalone_wrong_convention_assertion_refused(monkeypatch):
    monkeypatch.setattr(rmg_data_module, "database", None)
    with pytest.raises(PlasmaStateError, match="explicitly assert the ion thermo convention"):
        _initialize_proton(_load("electrochemical.py"), assertions={"H+": "electrochemical"})


def test_standalone_electron_remains_exempt(monkeypatch):
    monkeypatch.setattr(rmg_data_module, "database", None)
    library = _load("electrochemical.py")
    neutral = _species(library, "hydrogen", "H")
    electron = _electron()
    reactor = PlasmaReactor((300, "K"), (1, "bar"), {neutral: 1.0, electron: 1e-12},
                            (300, "K"), thermo_source_assertions=["e-"])
    reactor.initialize_model([neutral, electron], [], [], [])


def test_legacy_file_loader_resets_previous_convention():
    library = _load("ion.py")
    library.load_old(str(FIXTURES / "Dictionary.txt"), "", str(FIXTURES / "Library.txt"),
                     num_parameters=12, pattern=False)
    assert library.thermo_convention is None


def test_failed_file_load_cannot_retain_previous_convention(tmp_path):
    library = _load("ion.py")
    with pytest.raises(FileNotFoundError):
        library.load(str(tmp_path / "absent.py"))
    assert library.thermo_convention is None
