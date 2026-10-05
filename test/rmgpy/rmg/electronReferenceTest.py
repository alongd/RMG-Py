#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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

"""Complete material-reference census and pre-admission exception safety."""

import importlib.util
import inspect
import pickle
from pathlib import Path
from collections import namedtuple
from types import SimpleNamespace

import pytest

import rmgpy.kinetics as kinetics_module
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.library import KineticsLibrary, LibraryReaction
from rmgpy.exceptions import DatabaseError, NetworkError
from rmgpy.kinetics import (
    Arrhenius, EEDFChannel, Lindemann, MultiArrhenius, MultiPDepArrhenius,
    PDepArrhenius, Troe, ThirdBody,
)
from rmgpy.reaction import Reaction
from rmgpy.kinetics.model import KineticsModel
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.rmg.pdep import PDepNetwork
from rmgpy.species import Species
# Load the sibling helper by path, as plasmaKineticsConsolidationTest loads its helper.
_helpers_spec = importlib.util.spec_from_file_location(
    'electronTestHelpers', Path(__file__).with_name('electronTestHelpers.py'),
)
_helpers = importlib.util.module_from_spec(_helpers_spec)
_helpers_spec.loader.exec_module(_helpers)
(database, ORDERS, LIBRARIES, controlled_thermo, admission_state,
 registry_contents, _COLLIDER_BOUNDARIES) = (
    _helpers.database, _helpers.ORDERS, _helpers.LIBRARIES, _helpers.controlled_thermo,
    _helpers.admission_state, _helpers.registry_contents, _helpers._COLLIDER_BOUNDARIES,
)


KINETICS_CLASSES = tuple(sorted({cls for cls in vars(kinetics_module).values()
                                if inspect.isclass(cls) and issubclass(cls, KineticsModel)},
                               key=lambda cls: cls.__name__))
Cell = namedtuple('Cell', 'position holder child field form')


def make_rate(cls):
    if cls.__name__ == 'BadnellRRArrhenius':
        rate = cls(A=(1, 'm^3/(mol*s)'))
    elif cls.__name__ == 'VoronovEIArrhenius':
        rate = cls(A=(1, 'm^3/(mol*s)'), dE=1)
    elif cls.__name__ == 'Chebyshev':
        rate = cls(coeffs=[[1]], kunits='s^-1', Tmin=(300, 'K'), Tmax=(1000, 'K'),
                   Pmin=(0.1, 'bar'), Pmax=(10, 'bar'))
    elif cls.__name__ == 'KineticsData':
        rate = cls(Tdata=([300, 600], 'K'), kdata=([1, 2], 's^-1'))
    elif cls.__name__ == 'PDepKineticsData':
        rate = cls(Tdata=([300, 600], 'K'), Pdata=([1e5, 2e5], 'Pa'),
                   kdata=([[1, 2], [2, 3]], 's^-1'))
    else:
        rate = cls()
    if hasattr(rate, 'A') and rate.A is None:
        rate.A = (1, '') if cls.__name__.startswith('Sticking') else (1, 's^-1')
    for field in ('Ea', 'E0'):
        if hasattr(rate, field) and getattr(rate, field) is None:
            setattr(rate, field, (0, 'J/mol'))
    if hasattr(rate, 'arrhenius'):
        if cls.__name__ == 'MultiPDepArrhenius':
            member = make_rate(kinetics_module.PDepArrhenius)
        else:
            member = Arrhenius(A=(1, 's^-1'), n=0, Ea=(0, 'J/mol'))
        rate.arrhenius = [member]
        if hasattr(rate, 'pressures'):
            rate.pressures = ([1e5], 'Pa')
    if hasattr(rate, 'arrheniusLow'):
        rate.arrheniusLow = Arrhenius(A=(1, 'm^3/(mol*s)'), n=0, Ea=(0, 'J/mol'))
    if hasattr(rate, 'arrheniusHigh'):
        rate.arrheniusHigh = Arrhenius(A=(1, 's^-1'), n=0, Ea=(0, 'J/mol'))
    return rate


EFFICIENCY_CLASSES = tuple(cls for cls in KINETICS_CLASSES if hasattr(make_rate(cls), 'efficiencies'))
COVERAGE_CLASSES = tuple(cls for cls in KINETICS_CLASSES if hasattr(make_rate(cls), 'coverage_dependence'))
BEARING = [(cls, 'efficiencies', 'Molecule') for cls in EFFICIENCY_CLASSES]
BEARING += [(cls, 'coverage_dependence', form) for cls in COVERAGE_CLASSES
            for form in ('Species', 'Molecule')]


def census_cells():
    cells = []
    for cls in KINETICS_CLASSES:
        for field in ('reactants', 'products', 'pairs'):
            for form in ('Species', 'Molecule'):
                cells.append(Cell('reaction.' + field, cls, None, field, form))
        cells.append(Cell('reaction.specific_collider', cls, None, 'specific_collider', 'Species'))
    for cls, field, form in BEARING:
        cells.append(Cell('kinetics.' + field + '.key', cls, None, field, form))
    for parent in (kinetics_module.MultiArrhenius, kinetics_module.PDepArrhenius,
                   kinetics_module.MultiPDepArrhenius):
        for child, field, form in BEARING:
            cells.append(Cell('kinetics.arrhenius[].' + field + '.key', parent, child, field, form))
    for parent in EFFICIENCY_CLASSES:
        for child, field, form in BEARING:
            cells.append(Cell('kinetics.highPlimit.' + field + '.key', parent, child, field, form))
    for parent in (ThirdBody, Lindemann, Troe):
        for slot in ('arrheniusLow', 'arrheniusHigh'):
            if hasattr(make_rate(parent), slot):
                for form in ('Species', 'Molecule'):
                    cells.append(Cell('kinetics.' + slot + '.coverage_dependence.key', parent,
                                      kinetics_module.SurfaceArrhenius, slot, form))
    for form in ('Species', 'Molecule'):
        cells.append(Cell('reaction.network_kinetics.coverage_dependence.key', Arrhenius,
                          kinetics_module.SurfaceArrhenius, 'network_kinetics', form))
    for slot in ('SurfaceArrhenius', 'SurfaceChargeTransfer'):
        for child, field, form in BEARING:
            cells.append(Cell('reaction.' + slot + '.' + field + '.key', Arrhenius, child, field, form))
    return tuple(cells)


CENSUS_CELLS = census_cells()


def cell_id(cell):
    return '-'.join((cell.position, cell.holder.__name__,
                     cell.child.__name__ if cell.child else 'direct', cell.form))


def set_material_key(rate, field, reference):
    if field == 'efficiencies':
        rate.efficiencies = {reference: 2.0}
    else:
        rate.coverage_dependence = {reference: {'a': 0, 'm': 0, 'E': (0, 'J/mol')}}


def cell_reaction(cell, electron):
    source = Species(label='O2').from_smiles('[O][O]')
    product = Species(label='O').from_smiles('[O]')
    reference = (Species(label='electron_alias').from_adjacency_list('1 e u0 p0 c-1')
                 if electron else Species(label='ordinary_reference').from_smiles('N#N'))
    for species in (source, product, reference):
        controlled_thermo(species)
    held = reference if cell.form == 'Species' else reference.molecule[0]
    rate = make_rate(cell.holder)
    reaction = LibraryReaction(reactants=[source], products=[product, product],
                               kinetics=rate, reversible=False, library=LIBRARIES[0],
                               elementary_high_p=True)
    if cell.position.startswith('reaction.') and cell.child is None:
        if cell.field == 'pairs':
            reaction.pairs = [(held, product)]
        elif cell.field == 'specific_collider':
            reaction.specific_collider = held
        elif cell.field == 'reactants':
            reaction.reactants = [held]
        else:
            reaction.products = [held, product]
    elif cell.child is None:
        set_material_key(rate, cell.field, held)
    else:
        child = make_rate(cell.child)
        key_field = ('coverage_dependence' if cell.field in
                     ('arrheniusLow', 'arrheniusHigh', 'network_kinetics') else cell.field)
        set_material_key(child, key_field, held)
        if cell.position.startswith('kinetics.arrhenius[]'):
            rate.arrhenius = [child]
        elif cell.position.startswith('kinetics.highPlimit'):
            rate.highPlimit = child
        elif cell.position.startswith('reaction.'):
            slot = cell.position.split('.')[1]
            setattr(reaction, slot, child)
        else:
            setattr(rate, cell.field, child)
    return reaction


def admit_cell(model, reaction, cell):
    # Molecule participant normalization belongs to make_new_reaction; direct
    # registration accepts Species participant keys, and does not evaluate rates.
    if cell.field in ('reactants', 'products') and cell.form == 'Molecule' and cell.child is None:
        model.make_new_reaction(reaction, check_existing=False, generate_thermo=False,
                                generate_kinetics=False, perform_cut=False)
    else:
        model.register_reaction(reaction)


def ordinary_cell_snapshot(cell):
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    reaction = cell_reaction(cell, False)
    rate_before = pickle.dumps(reaction.kinetics)
    admit_cell(model, reaction, cell)
    assert pickle.dumps(reaction.kinetics) == rate_before
    return {'cell': cell_id(cell), 'registries': registry_contents(model)}


@pytest.mark.parametrize('wrapper', ('direct', 'multi', 'pdep', 'multi_pdep'))
def test_implicit_electron_eedf_channel_refuses_before_network_mutation(wrapper):
    """The persisted Ar => Ar* marker has no explicit electron participant."""
    source = Species(label='Ar').from_smiles('[Ar]')
    excited = Species(label='Ar*').from_smiles('[Ar]')
    channel = EEDFChannel('Ar -> Ar*', 'argon-lxcat-v1', 'ine')
    if wrapper == 'multi':
        rate = MultiArrhenius(arrhenius=[channel])
    elif wrapper == 'pdep':
        rate = PDepArrhenius(pressures=([1.0], 'bar'), arrhenius=[channel])
    elif wrapper == 'multi_pdep':
        rate = MultiPDepArrhenius(arrhenius=[
            PDepArrhenius(pressures=([1.0], 'bar'), arrhenius=[channel])])
    else:
        rate = channel
    reaction = Reaction(
        reactants=[source],
        products=[excited],
        kinetics=rate,
        electrons=0,
        reversible=False,
    )
    network = PDepNetwork(source=[source])
    before = dict(network.__dict__)

    with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks'):
        network.add_path_reaction(reaction)

    assert network.__dict__ == before


@pytest.mark.parametrize('order', [LIBRARIES])
@pytest.mark.parametrize('cell', CENSUS_CELLS, ids=cell_id)
@pytest.mark.parametrize('electron', (True, False), ids=('electron', 'ordinary'))
def test_species_reference_census_cell(database, order, cell, electron):
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    reaction = cell_reaction(cell, electron)
    if electron or cell.holder is EEDFChannel:
        before, candidate = admission_state(model), pickle.dumps(reaction)
        with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks') as error:
            admit_cell(model, reaction, cell)
        assert reaction.library in str(error.value)
        assert admission_state(model) == before and pickle.dumps(reaction) == candidate
    else:
        state = ordinary_cell_snapshot(cell)
        assert state['registries'][6]  # A reaction was registered without rate evaluation.


def efficiency_reaction(kind):
    reaction = cell_reaction(Cell('kinetics.efficiencies.key', kind, None, 'efficiencies', 'Molecule'), True)
    assert reaction.electrons == 0 and reaction.specific_collider is None
    assert not any(s.is_electron() for s in reaction.reactants + reaction.products)
    return reaction


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('kind', (ThirdBody, Lindemann, Troe), ids=lambda cls: cls.__name__)
@pytest.mark.parametrize('operation', _COLLIDER_BOUNDARIES +
                         ('seed', 'library', 'saved_seed', 'saved_library', 'pickle'))
def test_efficiency_electron_boundary_before_registration(
        database, order, kind, operation, monkeypatch, tmp_path):
    from arkane.pdep import PressureDependenceJob
    from rmgpy.pdep import Configuration
    import rmgpy.rmg.model as model_module
    import rmgpy.rmg.pdep as pdep_module

    monkeypatch.setattr(model_module, 'submit', controlled_thermo)
    reaction = efficiency_reaction(kind)
    library_name = reaction.library
    if operation in ('seed', 'library', 'saved_seed', 'saved_library'):
        library_name = 'electron-efficiency'
        stored = KineticsLibrary(label=library_name, name=library_name)
        stored.entries = {1: Entry(index=1, label=reaction.to_labeled_str(), item=reaction, data=reaction.kinetics)}
        if operation.startswith('saved_'):
            path = tmp_path / library_name
            path.mkdir()
            stored.save(str(path / 'reactions.py'))
            database.kinetics.load_libraries(str(tmp_path), libraries=[library_name], additive=True)
            reaction = database.kinetics.libraries[library_name].get_library_reactions()[0]
            assert reaction.elementary_high_p and reaction.electrons == 0
            assert any(key.is_electron() for key in reaction.kinetics.efficiencies)
        else:
            database.kinetics.libraries[library_name] = stored
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    source, isomer = reaction.reactants[0], reaction.products[0]
    network, other = PDepNetwork(source=[source]), PDepNetwork(source=[source])
    other.path_reactions = [reaction]
    job = PressureDependenceJob(network=None, Tmin=(300, 'K'), Tmax=(1000, 'K'), Tcount=2,
                               Pmin=(0.1, 'bar'), Pmax=(10, 'bar'), Pcount=2,
                               method='modified strong collision', interpolationModel=('Chebyshev', 2, 2))
    job.network, job.output_file = 'unchanged', None
    if operation in ('configurations', 'update', 'pickle'):
        network.path_reactions = [reaction]
        network.valid = True
    if operation == 'update':
        # A valid restored falloff path keeps its extracted high-pressure rate.
        reaction.network_kinetics = Arrhenius(A=(1, 's^-1'), n=0, Ea=(0, 'J/mol'))
    if operation == 'explore':
        # Exploration normalizes its isomer before producing the candidate.
        # Start with that normal graph state, as a library-loaded species does.
        for molecule in isomer.molecule:
            molecule.update()
        network.products = [Configuration(isomer)]
        monkeypatch.setattr(pdep_module, 'react_species', lambda reactants: [reaction])
    before, candidate = admission_state(model), pickle.dumps(reaction)
    net_before = {k: v[:] if isinstance(v, list) else v.copy() if isinstance(v, dict) else v
                  for k, v in network.__dict__.items()}
    with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks') as error:
        if operation == 'process':
            model.process_new_reactions([reaction], source, generate_thermo=False, generate_kinetics=False)
        elif operation == 'make_reaction':
            model.make_new_reaction(reaction, generate_thermo=False, generate_kinetics=False)
        elif operation == 'register':
            model.register_reaction(reaction)
        elif operation == 'make_pdep':
            model.make_new_pdep_reaction(reaction)
        elif operation == 'model_network':
            model.add_reaction_to_unimolecular_networks(reaction, source)
        elif operation == 'path':
            network.add_path_reaction(reaction)
        elif operation == 'merge':
            network.merge(other)
        elif operation == 'restore':
            network.__setstate__({'index': 999, 'path_reactions': [reaction]})
        elif operation == 'configurations':
            network.update_configurations(model)
        elif operation == 'update':
            network.update(model, job)
        elif operation == 'explore':
            network.explore_isomer(isomer)
        elif operation in ('library', 'saved_library'):
            model.add_reaction_library_to_edge(library_name)
        elif operation == 'pickle':
            pickle.loads(pickle.dumps(network))
        else:
            model.add_seed_mechanism_to_core(library_name)
    assert library_name in str(error.value)
    assert admission_state(model) == before and network.__dict__ == net_before
    assert job.network == 'unchanged' and pickle.dumps(reaction) == candidate


def source_failure_library(failure):
    from rmgpy.data.kinetics.common import format_external_library_provenance
    from rmgpy.reaction import Reaction

    species = Species(label='O2').from_smiles('[O][O]')
    reaction = Reaction(reactants=[species], products=[species], reversible=False)
    first = Arrhenius(A=(1, 's^-1'), comment='  first comment must survive  ')
    second = Arrhenius(A=(2, 's^-1'), comment='  second comment must survive  ')
    description = ('Originally from reaction library: A\n' + format_external_library_provenance('B', '/tmp/library-B')
                   if failure == 'DatabaseError' else
                   'Estimated using rate rule without brackets\nfamily: ordinary')
    stored = KineticsLibrary(label='source-conflict', auto_generated=True)
    stored.entries = {
        1: Entry(index=1, item=reaction, data=first,
                 long_desc='Estimated using rate rule [valid]\nfamily: ordinary'),
        2: Entry(index=2, item=reaction, data=second, long_desc=description),
    }
    return stored


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('operation', ('seed', 'library'))
@pytest.mark.parametrize('failure', ('DatabaseError', 'IndexError'))
def test_reconstruction_exception_restores_comments(database, order, operation, failure):
    stored = source_failure_library(failure)
    database.kinetics.libraries[stored.label] = stored
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    model.new_species_list, model.new_reaction_list = [object()], [object()]
    before = admission_state(model)
    comments = [entry.data.comment for entry in stored.entries.values()]
    method = model.add_seed_mechanism_to_core if operation == 'seed' else model.add_reaction_library_to_edge
    with pytest.raises(DatabaseError if failure == 'DatabaseError' else IndexError) as error:
        method(stored.label)
    expected = ('Conflicting external library provenance labels.' if failure == 'DatabaseError'
                else 'list index out of range')
    assert str(error.value) == expected
    assert [entry.data.comment for entry in stored.entries.values()] == comments
    assert admission_state(model) == before


@pytest.mark.parametrize('order', [LIBRARIES])
@pytest.mark.parametrize('electron', (True, False), ids=('electron', 'ordinary'))
@pytest.mark.parametrize('position', ('reaction', 'kinetics', 'wrapper', 'dictionary_key',
                                     'dictionary_value', 'slots', 'class_field', 'object_array',
                                     'provenance_wrapper'))
def test_unclassified_species_reference_fails_closed(database, order, electron, position):
    import numpy as np

    class ExtendedArrhenius(Arrhenius):
        pass

    class SlotCarrier:
        __slots__ = ('material',)

    reaction = cell_reaction(Cell('reaction.specific_collider', Arrhenius, None,
                                  'specific_collider', 'Species'), electron)
    reference, reaction.specific_collider = reaction.specific_collider, None
    reaction.kinetics = ExtendedArrhenius(A=(1, 's^-1'), n=0, Ea=(0, 'J/mol'))
    if position == 'reaction':
        reaction.unclassified_material = reference
    elif position == 'class_field':
        ExtendedArrhenius.unclassified_material = reference
    elif position == 'wrapper':
        reaction.kinetics.unclassified_material = SimpleNamespace(material=reference)
    elif position == 'dictionary_key':
        reaction.kinetics.unclassified_material = {reference: 2}
    elif position == 'dictionary_value':
        reaction.kinetics.unclassified_material = {'material': reference}
    elif position == 'slots':
        holder = SlotCarrier()
        holder.material = reference
        reaction.kinetics.unclassified_material = holder
    elif position == 'object_array':
        reaction.kinetics.unclassified_material = np.array([reference], dtype=object)
    elif position == 'provenance_wrapper':
        reaction.entry = SimpleNamespace(material=reference)
    else:
        reaction.kinetics.unclassified_material = reference
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    before = admission_state(model)
    with pytest.raises(NetworkError, match='Unclassified species reference') as error:
        model.register_reaction(reaction)
    assert reaction.library in str(error.value) and admission_state(model) == before


@pytest.mark.parametrize('order', [LIBRARIES])
@pytest.mark.parametrize('electron', (True, False), ids=('electron', 'ordinary'))
def test_generic_nested_kinetics_and_cycles(database, order, electron):
    from rmgpy.rmg.pdep import _has_electron_participant, _reaction_species_references

    class ExtendedArrhenius(Arrhenius):
        pass

    reaction = cell_reaction(Cell('kinetics.efficiencies.key', Troe, None,
                                  'efficiencies', 'Molecule'), electron)
    leaf = reaction.kinetics
    parent = ExtendedArrhenius(A=(1, 's^-1'), n=0, Ea=(0, 'J/mol'))
    holder = SimpleNamespace(rate=leaf)
    holder.parent = parent
    parent.new_member_slot = {'expressions': [holder]}
    parent.self_reference = parent
    reaction.kinetics = parent
    with pytest.raises(ValueError, match='exact type'):
        list(_reaction_species_references(reaction))
    assert _has_electron_participant(reaction)
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    before = admission_state(model)
    with pytest.raises(NetworkError):
        model.register_reaction(reaction)
    assert admission_state(model) == before


@pytest.mark.parametrize('order', [LIBRARIES])
@pytest.mark.parametrize('electron', (True, False), ids=('electron', 'ordinary'))
def test_unknown_class_refuses_before_later_fields(database, order, electron):
    from rmgpy.rmg.pdep import _has_electron_participant

    class ExtendedArrhenius(Arrhenius):
        pass

    reaction = cell_reaction(Cell('reaction.specific_collider', Arrhenius, None,
                                  'specific_collider', 'Species'), electron)
    reference, reaction.specific_collider = reaction.specific_collider, None
    reaction.kinetics = ExtendedArrhenius(A=(1, 's^-1'), n=0, Ea=(0, 'J/mol'))
    assert _has_electron_participant(reaction)
    ExtendedArrhenius.later_material = reference
    assert _has_electron_participant(reaction)
    del ExtendedArrhenius.later_material
    assert _has_electron_participant(reaction)
    reaction.kinetics.later_material = reference
    assert _has_electron_participant(reaction)


class ArbitraryPreAdmissionFailure(BaseException):
    pass


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('operation', ('seed', 'library'))
@pytest.mark.parametrize('phase', ('reconstruction', 'preflight'))
@pytest.mark.parametrize('exception_class', (RuntimeError, ArbitraryPreAdmissionFailure))
def test_arbitrary_exception_restores_comments(
        database, order, operation, phase, exception_class, monkeypatch):
    stored = source_failure_library('IndexError')
    stored.entries[2].long_desc = 'Estimated using rate rule [valid]\nfamily: ordinary'
    database.kinetics.libraries[stored.label] = stored
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    model.new_species_list, model.new_reaction_list = [object()], [object()]
    before = admission_state(model)
    comments = [entry.data.comment for entry in stored.entries.values()]
    sentinel = exception_class('original exception', object())

    def fail(*args, **kwargs):
        if phase == 'reconstruction':
            stored.entries[1].data.comment = comments[0].strip()
        raise sentinel

    if phase == 'reconstruction':
        monkeypatch.setattr(stored, 'get_library_reactions', fail)
    else:
        monkeypatch.setattr(model, '_preflight_electron_routing', fail)
    method = model.add_seed_mechanism_to_core if operation == 'seed' else model.add_reaction_library_to_edge
    with pytest.raises(exception_class) as error:
        method(stored.label)
    assert error.value is sentinel
    assert [entry.data.comment for entry in stored.entries.values()] == comments
    assert admission_state(model) == before


class ExtendedReferenceArrhenius(Arrhenius):
    """A Python extension with state outside the native kinetics schema."""


class ReferenceList(list):
    pass


class ReferenceDict(dict):
    pass


class ReferenceSlots(list):
    __slots__ = ('material',)


class CallableReference:
    __slots__ = ('material',)

    def __call__(self):
        return None


@pytest.mark.parametrize('order', [LIBRARIES])
@pytest.mark.parametrize('reference_kind', ('electron', 'ordinary', 'numeric'))
@pytest.mark.parametrize('carrier', (
    'list', 'dict', 'array', 'slots', 'callable_class', 'callable_function',
    'static_callable', 'class_callable', 'closure_callable',
))
@pytest.mark.parametrize('boundary', ('register', 'process', 'path'))
def test_container_and_callable_extension_state_fails_closed(
        database, order, reference_kind, carrier, boundary, monkeypatch):
    import numpy as np
    import rmgpy.rmg.model as model_module

    class ReferenceArray(np.ndarray):
        pass

    monkeypatch.setattr(model_module, 'submit', controlled_thermo)
    reaction = cell_reaction(CENSUS_CELLS[0], False)
    reaction.kinetics = ExtendedReferenceArrhenius(A=(1, 's^-1'), n=0, Ea=(0, 'J/mol'))
    reference = (42 if reference_kind == 'numeric' else
                 Species().from_adjacency_list('1 e u0 p0 c-1') if reference_kind == 'electron'
                 else Species().from_smiles('N#N'))

    def function_carrier(self):
        return None

    factories = {
        'callable_function': lambda: function_carrier, 'list': ReferenceList,
        'dict': ReferenceDict, 'slots': ReferenceSlots, 'callable_class': CallableReference,
        'static_callable': CallableReference, 'class_callable': CallableReference,
    }
    holder = factories.get(carrier, lambda: np.array([1.0]).view(ReferenceArray))()
    if carrier == 'closure_callable':
        holder = lambda self: reference
    else:
        holder.material = reference
    if carrier in ('static_callable', 'class_callable'):
        descriptor = staticmethod(holder) if carrier == 'static_callable' else classmethod(holder)
        monkeypatch.setattr(ExtendedReferenceArrhenius, 'material_carrier', descriptor, raising=False)
    elif carrier in ('callable_class', 'callable_function', 'closure_callable'):
        monkeypatch.setattr(ExtendedReferenceArrhenius, 'material_carrier', holder, raising=False)
    else:
        reaction.kinetics.material_carrier = holder
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    network = PDepNetwork(source=reaction.reactants)
    before = admission_state(model)
    with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks') as error:
        if boundary == 'register':
            model.register_reaction(reaction)
        elif boundary == 'process':
            model.process_new_reactions([reaction], reaction.reactants[0],
                                        generate_thermo=False, generate_kinetics=False)
        else:
            network.add_path_reaction(reaction)
    assert reaction.library in str(error.value)
    assert admission_state(model) == before
    assert network.path_reactions == [] and network.net_reactions == []


def restore_census_rate(rate, state):
    """Restore the material fields in a complete legacy cache fixture."""
    for name, value in state.items():
        setattr(rate, name, value)


class CensusCachePickler(pickle.Pickler):
    """Emit complete cached material state, without changing global reducers.

    Native reducers omit several artificial census slots (highPlimit and the
    never-assigned reaction type slots); PDepKineticsModel's native constructor
    tuple also does not round-trip. A cache fixture must retain these positions
    to exercise routing of a restored graph, rather than pass by losing its
    electron. The independent r144 regressions use native pickle unchanged.
    """

    def reducer_override(self, value):
        from rmgpy.rmg.pdep import PDepReaction

        if isinstance(value, KineticsModel):
            fields = ('efficiencies', 'coverage_dependence', 'highPlimit',
                      'arrhenius', 'arrheniusLow', 'arrheniusHigh')
            state = {name: getattr(value, name) for name in fields if hasattr(value, name)}
            return (make_rate, (type(value),), state, None, None, restore_census_rate)
        if isinstance(value, PDepReaction):
            reduction = value.__reduce__()
            state = dict(reduction[2])
            for name in ('SurfaceArrhenius', 'SurfaceChargeTransfer'):
                state[name] = getattr(value, name)
            return reduction[0], reduction[1], state
        return NotImplemented


def census_cache_bytes(root):
    import io

    buffer = io.BytesIO()
    CensusCachePickler(buffer).dump(root)
    return buffer.getvalue()


def circular_census_root(cell, electron, root_kind, slot):
    from rmgpy.data.kinetics.family import _NOT_REPRODUCED, apply_reaction_state, reaction_state, state_fields
    from rmgpy.rmg.model import ReactionModel
    from rmgpy.rmg.pdep import PDepReaction

    original = cell_reaction(cell, electron)
    reaction = PDepReaction()
    apply_reaction_state(reaction, reaction_state(original, _NOT_REPRODUCED, state_fields(original)))
    for name in ('SurfaceArrhenius', 'SurfaceChargeTransfer'):
        setattr(reaction, name, getattr(original, name))
    network = PDepNetwork(source=reaction.reactants)
    reaction.network = network
    setattr(network, slot, [reaction])  # A cached legacy graph, before admission.
    if root_kind == 'network':
        root = network
    elif root_kind == 'reaction':
        root = reaction
    elif root_kind == 'reaction_model':
        root = ReactionModel(reactions=[reaction])
    else:
        root = CoreEdgeReactionModel()
        root.core.reactions = [reaction]
        root.network_list = [network]
    return root


@pytest.mark.parametrize('cell', CENSUS_CELLS, ids=cell_id)
@pytest.mark.parametrize('slot', ('path_reactions', 'net_reactions'))
@pytest.mark.parametrize('root_kind', ('network', 'reaction', 'reaction_model', 'core_edge_model'))
@pytest.mark.parametrize('electron', (True, False), ids=('electron', 'ordinary'))
def test_circular_pickle_census_first_use(cell, slot, root_kind, electron):
    root = circular_census_root(cell, electron, root_kind, slot)
    data = census_cache_bytes(root)

    def restored_network():
        restored = pickle.loads(data)
        if root_kind == 'network':
            return restored
        reaction = (restored if root_kind == 'reaction' else restored.reactions[0]
                    if root_kind == 'reaction_model' else restored.core.reactions[0])
        return reaction.network

    if electron or cell.holder is EEDFChannel:
        with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks') as error:
            network = restored_network()
            # The complete graph is available here, including reactions whose
            # state was installed after the network's own __setstate__.
            network.add_path_reaction(getattr(network, slot)[0])
        assert LIBRARIES[0] in str(error.value)
    else:
        network = restored_network()
        restored_reaction = getattr(network, slot)[0]
        assert restored_reaction.network is network
        network.add_path_reaction(restored_reaction)


@pytest.mark.parametrize('slot', ('path_reactions', 'net_reactions'))
def test_restored_network_checks_before_inherited_numerical_mutation(slot):
    cell = Cell('reaction.reactants', Arrhenius, None, 'reactants', 'Species')
    reaction = pickle.loads(pickle.dumps(circular_census_root(cell, True, 'reaction', slot)))
    network = reaction.network
    before = dict(network.__dict__)
    with pytest.raises(NetworkError, match=LIBRARIES[0]):
        network.initialize(300, 1000, 1e4, 1e6)
    assert network.__dict__ == before


def test_restored_network_does_not_cache_placeholder_absence():
    from rmgpy.rmg.pdep import PDepReaction

    reaction = PDepReaction(reactants=[], products=[])
    reaction.library = LIBRARIES[0]
    network = PDepNetwork()
    network.__setstate__({'path_reactions': [reaction], 'net_reactions': []})
    assert network.path_reactions == [reaction]  # Inspection is allowed.
    network._check_reactions()  # A placeholder check cannot authorize later use.
    reaction.reactants = [Species().from_adjacency_list('1 e u0 p0 c-1')]
    with pytest.raises(NetworkError, match=LIBRARIES[0]):
        network.get_leak_coefficient(300, 1e5)


@pytest.mark.parametrize('form', ('reactants', 'efficiencies'))
@pytest.mark.parametrize('slot', ('path_reactions', 'net_reactions'))
@pytest.mark.parametrize('root_kind', ('network', 'reaction', 'reaction_model', 'core_edge_model'))
@pytest.mark.parametrize('electron', (True, False), ids=('electron', 'ordinary'))
def test_circular_native_pickle_review_witness(form, slot, root_kind, electron):
    cell = (Cell('reaction.reactants', Arrhenius, None, 'reactants', 'Species')
            if form == 'reactants' else
            Cell('kinetics.efficiencies.key', Lindemann, None, 'efficiencies', 'Molecule'))
    data = pickle.dumps(circular_census_root(cell, electron, root_kind, slot))

    def first_use():
        restored = pickle.loads(data)
        if root_kind == 'network':
            network = restored
        else:
            reaction = (restored if root_kind == 'reaction' else restored.reactions[0]
                        if root_kind == 'reaction_model' else restored.core.reactions[0])
            network = reaction.network
        reaction = getattr(network, slot)[0]
        network.add_path_reaction(reaction)
        return network, reaction

    if electron:
        with pytest.raises(NetworkError, match=LIBRARIES[0]):
            first_use()
    else:
        network, reaction = first_use()
        assert reaction.network is network


@pytest.mark.parametrize('electron', (True, False), ids=('electron', 'ordinary'))
@pytest.mark.parametrize('position', ('reactants', 'coverage_dependence'))
def test_species_material_representations_are_all_inspected(electron, position):
    cell = (Cell('reaction.reactants', Arrhenius, None, 'reactants', 'Species')
            if position == 'reactants' else
            Cell('kinetics.coverage_dependence.key', kinetics_module.SurfaceArrhenius,
                 None, 'coverage_dependence', 'Species'))
    reaction = cell_reaction(cell, False)
    reference = (reaction.reactants[0] if position == 'reactants'
                 else next(iter(reaction.kinetics.coverage_dependence)))
    other = (Species().from_adjacency_list('1 e u0 p0 c-1') if electron
             else Species().from_smiles('N#N'))
    reference.molecule.append(other.molecule[0])
    assert not reference.is_electron()  # The native predicate looks at its first molecule.
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    if electron:
        before = admission_state(model)
        with pytest.raises(NetworkError, match=LIBRARIES[0]):
            model.register_reaction(reaction)
        assert admission_state(model) == before
    else:
        model.register_reaction(reaction)
        assert model.reaction_dict


@pytest.mark.parametrize('electron', (True, False), ids=('electron', 'ordinary'))
def test_method_closure_species_data_fails_closed(electron):
    reference = (Species().from_adjacency_list('1 e u0 p0 c-1') if electron
                 else Species().from_smiles('N#N'))

    class CapturingRate(Arrhenius):
        def stored_reference(self):
            return reference

    reaction = cell_reaction(CENSUS_CELLS[0], False)
    reaction.kinetics = CapturingRate(A=(1, 's^-1'), n=0, Ea=(0, 'J/mol'))
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    before = admission_state(model)
    with pytest.raises(NetworkError, match=LIBRARIES[0]):
        model.register_reaction(reaction)
    assert admission_state(model) == before


def test_numeric_method_closure_subclass_refuses():
    coefficient = 42

    class OrdinaryRate(Arrhenius):
        def ordinary_method(self):
            assert super().is_temperature_valid(300)
            return coefficient

    reaction = cell_reaction(CENSUS_CELLS[0], False)
    reaction.kinetics = OrdinaryRate(A=(1, 's^-1'), n=0, Ea=(0, 'J/mol'))
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    with pytest.raises(NetworkError, match='exact type'):
        model.register_reaction(reaction)
    assert not model.reaction_dict and reaction.kinetics.ordinary_method() == 42
