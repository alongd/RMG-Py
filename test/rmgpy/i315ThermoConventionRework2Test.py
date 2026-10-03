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

"""Regression tests for complete edge validation and isolated library parsing."""

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
from rmgpy.thermo import ThermoData


DATA = Path(__file__).parent / 'test_data'
FIXTURES = DATA / 'i315_thermo_convention'
NEW_FIXTURES = DATA / 'i315_thermo_convention_rework2'


def _load(path, context=None):
    context = ThermoDatabase() if context is None else context
    return ThermoLibrary(label=path.stem).load(
        str(path), context.local_context, context.global_context)


def _species(library, entry, label):
    source = library.entries[entry]
    return Species(label=label, molecule=[copy.deepcopy(source.item)],
                   thermo=copy.deepcopy(source.data))


def _database(monkeypatch, library):
    thermo = ThermoDatabase()
    thermo.library_order = [library.label]
    thermo.libraries = {library.label: library}
    monkeypatch.setattr(rmg_data_module, 'database', SimpleNamespace(thermo=thermo))


def _edge_fixture():
    actual = _load(FIXTURES / 'electrochemical.py')
    proton = _species(actual, 'proton', 'H+')
    oxide = _species(actual, 'oxide', 'O-')
    hydrogen = _species(actual, 'hydrogen', 'H')
    oxygen = _species(actual, 'oxygen', 'O')
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    core = [hydrogen, oxygen, electron]
    reactor = PlasmaReactor((300, 'K'), (1, 'bar'),
                            {hydrogen: 0.5, oxygen: 0.5, electron: 1e-12}, (300, 'K'))
    reaction = Reaction(reactants=[proton, oxide], products=[hydrogen, oxygen],
                        reversible=True, kinetics=Arrhenius(A=(1, 'm^3/(mol*s)')))
    return reactor, core, proton, oxide, reaction


@pytest.mark.parametrize('reference', ['ion.py', 'reviewed_legacy.py'])
def test_unmatched_charged_edge_refused_before_reverse_growth_rates(monkeypatch, reference):
    """An ion-declared mismatch and no candidate must both refuse the bad proton."""
    library = _load(FIXTURES / reference)
    assert library.thermo_convention == 'ion'
    _database(monkeypatch, library)
    reactor, core, proton, oxide, reaction = _edge_fixture()
    with pytest.raises(PlasmaStateError, match=r'edge species.*H\+.*could not be value-matched'):
        reactor.initialize_model(core, [], [proton, oxide], [reaction])
    assert reactor.num_edge_reactions == -1
    assert not reactor._plasma_validated


def test_missing_charged_edge_thermo_is_refused(monkeypatch):
    _database(monkeypatch, _load(FIXTURES / 'ion.py'))
    reactor, core, proton, _, _ = _edge_fixture()
    proton.thermo = None
    with pytest.raises(PlasmaStateError, match=r'edge species.*H\+.*has no thermo data'):
        reactor.initialize_model(core, [], [proton], [])


def test_malformed_ion_reference_does_not_admit_edge_thermo(monkeypatch):
    _database(monkeypatch, _load(NEW_FIXTURES / 'malformed_ion.py'))
    reactor, core, proton, _, _ = _edge_fixture()
    with pytest.raises(PlasmaStateError, match=r'edge species.*H\+.*malformed'):
        reactor.initialize_model(core, [], [proton], [])


def test_valid_ion_edge_and_unmatched_neutral_edge_are_accepted(monkeypatch):
    library = _load(FIXTURES / 'ion.py')
    _database(monkeypatch, library)
    reactor, core, _, _, _ = _edge_fixture()
    proton = _species(library, 'proton', 'H+')
    neutral = Species(label='estimated-neutral', smiles='N#N')
    neutral.thermo = copy.deepcopy(core[0].thermo)
    reactor.initialize_model(core, [], [proton, neutral], [])
    assert reactor.thermo_provenance_diagnostics['H+'].startswith('value-matched to library ion/proton')
    assert 'estimated-neutral' not in reactor.thermo_provenance_diagnostics


def test_earlier_library_cannot_corrupt_pinned_thermo_values():
    reviewed = FIXTURES / 'reviewed_legacy.py'
    clean = _load(reviewed)
    expected = clean.entries['[Arp]'].data.H298.value_si
    context = ThermoDatabase()
    context.load_libraries(str(NEW_FIXTURES), libraries=[
        str(NEW_FIXTURES / 'namespace_poison.py'), str(reviewed)])
    loaded = context.libraries['reviewed_legacy']
    assert loaded.thermo_convention == 'ion'
    assert loaded._loaded_file_sha256 == clean._loaded_file_sha256
    assert loaded.entries['[Arp]'].data.H298.value_si == pytest.approx(expected)
    assert expected > 1.5e6
    assert context.local_context['ThermoData'] is ThermoData
    assert 'parser_leak' not in context.global_context


def test_library_execution_does_not_mutate_caller_parser_dictionaries():
    context = ThermoDatabase()
    context.local_context['user_value'] = 'preserve this caller binding'
    context.global_context['global_value'] = 'preserve this caller binding'
    locals_before = dict(context.local_context)
    globals_before = dict(context.global_context)
    _load(NEW_FIXTURES / 'namespace_poison.py', context)
    assert context.local_context == locals_before
    assert context.global_context == globals_before


def _corrupt_constructor(*args, **kwargs):
    return ThermoData(Tdata=([300, 1500], 'K'), Cpdata=([20, 20], 'J/(mol*K)'),
                      H298=(3, 'J/mol'), S298=(0, 'J/(mol*K)'))


def test_pinned_library_uses_trusted_constructors_even_with_poisoned_caller_context():
    reviewed = FIXTURES / 'reviewed_legacy.py'
    clean = _load(reviewed)
    context = ThermoDatabase()
    for name in ['ThermoData', 'Wilhoit', 'NASAPolynomial', 'NASA']:
        context.local_context[name] = _corrupt_constructor
        context.global_context[name] = _corrupt_constructor
    loaded = _load(reviewed, context)
    for label in clean.entries:
        assert type(loaded.entries[label].data) is type(clean.entries[label].data)
        assert loaded.entries[label].data.get_enthalpy(298.15) == pytest.approx(
            clean.entries[label].data.get_enthalpy(298.15))
    assert context.local_context['ThermoData'] is _corrupt_constructor


def test_refused_reinitialization_clears_rate_validation_latch(monkeypatch):
    _database(monkeypatch, _load(FIXTURES / 'ion.py'))
    reactor, core, proton, _, _ = _edge_fixture()
    reactor.initialize_model(core, [], [], [])
    assert reactor._plasma_validated
    # The existing core guard already refuses this bad proton at the base.
    with pytest.raises(PlasmaStateError, match='could not be value-matched'):
        reactor.initialize_model(core + [proton], [], [], [])
    assert not reactor._plasma_validated
    with pytest.raises(PlasmaStateError, match='before.*validated'):
        reactor.generate_rate_coefficients([], [])


class _LateRefusalReactor(PlasmaReactor):
    def _resolve_wall_state(self, core_species):
        raise PlasmaStateError('test refusal after the rate-validation latch was set')


def test_late_initialization_refusal_also_clears_rate_validation_latch(monkeypatch):
    monkeypatch.setattr(rmg_data_module, 'database', None)
    neutral = Species(label='H').from_adjacency_list('multiplicity 2\n1 H u1 p0 c0')
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    reactor = _LateRefusalReactor((300, 'K'), (1, 'bar'),
                                  {neutral: 1.0, electron: 1e-12}, (300, 'K'))
    with pytest.raises(PlasmaStateError, match='test refusal'):
        reactor.initialize_model([neutral, electron], [], [], [])
    assert not reactor._plasma_validated
    with pytest.raises(PlasmaStateError, match='before.*validated'):
        reactor.generate_rate_coefficients([], [])


@pytest.mark.parametrize('job', ['plasma-only', 'mixed', 'seeded-lithium', 'declared-zero-lithium'])
def test_impossible_library_element_excluded_before_edge_thermo(monkeypatch, caplog, job):
    """An unrelated lithium library entry cannot populate an argon-only plasma edge."""
    import logging
    import rmgpy.rmg.input as input_module
    from rmgpy.rmg.model import CoreEdgeReactionModel

    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    argon = Species(label='Ar', smiles='[Ar]')
    argon_ion = Species(label='Arp').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    lithium = Species(label='Li', smiles='[Li]')
    lithium_ion = Species(label='Lip').from_adjacency_list('1 Li u0 p0 c+1')
    argon_reaction = Reaction(reactants=[argon_ion, electron], products=[argon], reversible=False)
    lithium_reaction = Reaction(reactants=[lithium_ion, electron], products=[lithium], reversible=False)
    library = SimpleNamespace(name='mixed-elements', get_library_reactions=lambda: [lithium_reaction, argon_reaction])
    kinetics = SimpleNamespace(families={}, resolve_library=lambda _: library)
    monkeypatch.setattr(rmg_data_module, 'database', SimpleNamespace(kinetics=kinetics))
    initial_species = [argon, electron]
    if job == 'declared-zero-lithium':
        initial_species.append(lithium)
    monkeypatch.setattr(input_module, 'rmg', SimpleNamespace(
        initial_species=initial_species, job_is_all_plasma=lambda: job != 'mixed'))
    model = CoreEdgeReactionModel()
    model.save_edge_species = False
    if job == 'seeded-lithium':
        model.core.species.append(lithium)
    imported = []

    def record_before_species_or_thermo_creation(reaction):
        imported.append(reaction)
        return reaction, False

    monkeypatch.setattr(model, 'make_new_reaction', record_before_species_or_thermo_creation)
    with caplog.at_level(logging.INFO):
        model.add_reaction_library_to_edge('mixed-elements')
    if job == 'plasma-only':
        assert imported == [argon_reaction]
        assert 'Excluding reaction' in caplog.text
        assert 'mixed-elements' in caplog.text and 'Li' in caplog.text
        assert 'absent from the input and seed inventory' in caplog.text
    else:
        assert imported == [lithium_reaction, argon_reaction]
        assert 'Excluding reaction' not in caplog.text
